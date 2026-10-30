"""Resumable TerrainRGB tiles and verified PMTiles packaging from numeric DEMs."""
import hashlib
import io
import json
import math
import multiprocessing
import os
import sqlite3
import tempfile
import time
from contextlib import ExitStack, closing
from concurrent.futures import FIRST_COMPLETED,ProcessPoolExecutor,wait
from pathlib import Path

import numpy as np
import rasterio
from PIL import Image
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from rasterio.warp import reproject,transform_bounds
from rasterio.windows import Window,transform as window_transform

from terrain_rgb import encode_terrain_rgb

HALF_WORLD = math.pi * 6378137
TILE_SIZE = 512
QUANTIZATION = '16-max(z,5)'
ROWS_PER_TASK = 32


def round_bits(zoom):
    if not 0 <= zoom <= 12:
        raise ValueError('Production tiles support zooms 0 through 12')
    return 16-max(zoom,5)


def grid(bounds, zoom, tile_size=TILE_SIZE):
    west,south,east,north = transform_bounds('EPSG:4326','EPSG:3857',*bounds,densify_pts=21)
    span = 2 * HALF_WORLD / 2**zoom
    xmin=max(0,math.floor((west+HALF_WORLD)/span))
    xmax=min(2**zoom-1,math.ceil((east+HALF_WORLD)/span)-1)
    ymin=max(0,math.floor((HALF_WORLD-north)/span))
    ymax=min(2**zoom-1,math.ceil((HALF_WORLD-south)/span)-1)
    if not -HALF_WORLD <= south < north <= HALF_WORLD:
        raise ValueError('Bounds exceed Web Mercator')
    return dict(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
                transform=from_origin(-HALF_WORLD+xmin*span,HALF_WORLD-ymin*span,
                                      span/tile_size,span/tile_size),
                width=(xmax-xmin+1)*tile_size,height=(ymax-ymin+1)*tile_size,
                resolution=span/tile_size)


class Sampler:
    """Memory-bounded numeric and validity warps on one shared output grid."""
    def __init__(self, heights, validity, zoom, tile_size=TILE_SIZE, bounds=None):
        self.stack=ExitStack()
        self.tile_size=tile_size
        self.zoom=zoom
        self.stack.enter_context(rasterio.Env(GDAL_CACHEMAX=64*1024**2,GDAL_NUM_THREADS='1',
                                             VRT_NUM_THREADS='1',GDAL_MAX_DATASET_POOL_SIZE=48))
        source=self.source=self.stack.enter_context(rasterio.open(heights))
        mask=self.validity=self.stack.enter_context(rasterio.open(validity))
        if source.nodata is None or source.count!=1 or source.dtypes!=('float32',):
            raise ValueError('Expected Float32 elevations with declared nodata')
        if (source.transform!=mask.transform or source.shape!=mask.shape
                or source.crs!=mask.crs or mask.nodata is not None):
            raise ValueError('Validity must share the source grid and treat zero as data')
        self.grid=grid(bounds or source.bounds,zoom,tile_size)
        self.method=Resampling.average if zoom<=11 else Resampling.bilinear

    def read(self,x,y):
        w=Window((x-self.grid['xmin'])*self.tile_size,(y-self.grid['ymin'])*self.tile_size,
                 self.tile_size,self.tile_size)
        common=dict(src_crs='EPSG:4326',dst_crs='EPSG:3857',
                    dst_transform=window_transform(w,self.grid['transform']),
                    warp_mem_limit=32,num_threads=1,tolerance=1e-6)
        # WarpedVRT block reads can allocate the entire native source window
        # at low zoom, regardless of warp_mem_limit. reproject uses GDAL's
        # chunking warp operation, which subdivides that window to fit.
        mask=np.zeros((self.tile_size,self.tile_size),dtype=np.uint8)
        reproject(rasterio.band(self.validity,1),mask,**common,
                  src_nodata=None,dst_nodata=0,resampling=Resampling.nearest)
        valid=mask==255
        if not valid.any():
            return np.full(valid.shape,np.nan,dtype=np.float32),valid
        h=np.full(valid.shape,np.nan,dtype=np.float32)
        if self.zoom<5:
            # A world-sized destination includes the antimeridian even though
            # this DEM is regional. GDAL can request a full-width native source
            # window there, exceeding its warp budget at the minimum chunk size.
            # Only valid target pixels will be encoded: crop the numeric warp
            # to their bounding rectangle, preserving the exact output grid.
            rows,cols=np.nonzero(valid)
            top,bottom=int(rows.min()),int(rows.max())+1
            left,right=int(cols.min()),int(cols.max())+1
            cropped=np.full((bottom-top,right-left),np.nan,dtype=np.float32)
            cropped_common={**common,'dst_transform':window_transform(
                Window(left,top,right-left,bottom-top),common['dst_transform'])}
            reproject(rasterio.band(self.source,1),cropped,**cropped_common,
                      src_nodata=self.source.nodata,dst_nodata=np.nan,resampling=self.method)
            h[top:bottom,left:right]=cropped
        else:
            reproject(rasterio.band(self.source,1),h,**common,
                      src_nodata=self.source.nodata,dst_nodata=np.nan,resampling=self.method)
        valid &= np.isfinite(h)
        return h,valid

    def __enter__(self):
        return self

    def __exit__(self,*args):
        self.stack.close()


def webp_bytes(rgba):
    stream=io.BytesIO()
    Image.fromarray(rgba).save(stream,format='WEBP',lossless=True,exact=True,method=4)
    data=stream.getvalue()
    with Image.open(io.BytesIO(data)) as decoded:
        if not np.array_equal(np.asarray(decoded.convert('RGBA')),rgba):
            raise ValueError('WebP changed RGBA bytes')
    return data


def validate_rgba(rgba,zoom,valid):
    """Validate display bytes against the explicit input mask, never inferred RGB."""
    alpha=rgba[...,3]
    if not np.all(alpha==255):
        raise ValueError('Unexpected alpha value')
    valid=np.asarray(valid,dtype=bool)
    if valid.shape!=alpha.shape:
        raise ValueError('Validity mask must match the tile shape')
    if np.any(rgba[~valid]!=np.array([1,134,160,255],dtype=np.uint8)):
        raise ValueError('Nodata RGBA placeholder changed')
    rgb=rgba[valid,:3].astype(np.int32)
    codes=rgb[:,0]*65536+rgb[:,1]*256+rgb[:,2]
    if np.any(codes % 2**round_bits(zoom)):
        raise ValueError('Valid heights do not obey quantization schedule')
    return int(valid.sum())


def encode_rows(sampler,empty,z,x,start,end):
    """One bounded task; workers never access the destination SQLite database."""
    from pmtiles.tile import zxy_to_tileid
    rows=[]
    for y in range(start,end):
        h,valid=sampler.read(x,y);count=int(valid.sum());error=0.0
        if count:
            rgba=encode_terrain_rgb(h,valid=valid,round_bits=round_bits(z))
            validate_rgba(rgba,z,valid)
            rgb=rgba[valid,:3].astype(np.int32)
            decoded=-10000+.1*(rgb[:,0]*65536+rgb[:,1]*256+rgb[:,2])
            error=float(np.abs(decoded-h[valid]).max())
            if error>.05*2**round_bits(z)+1e-5:
                raise ValueError('Encoding exceeded quantization tolerance')
            data=webp_bytes(rgba)
        else:data=empty
        rows.append((zxy_to_tileid(z,x,y),x,y,data,hashlib.sha256(data).hexdigest(),count,error))
    return x,rows


def init_tile_worker(heights,validity,z,bounds,empty):
    import atexit
    global worker_sampler,worker_empty,worker_zoom
    worker_sampler=Sampler(heights,validity,z,bounds=bounds)
    worker_empty=empty;worker_zoom=z
    atexit.register(worker_sampler.stack.close)


def run_tile_task(task):
    x,start,end=task
    return encode_rows(worker_sampler,worker_empty,worker_zoom,x,start,end)


def parallel_rows(heights,validity,z,bounds,empty,tasks,workers):
    # Spawn avoids inheriting open GDAL datasets. Each worker owns its readers,
    # 64 MiB block cache and one 32 MiB, single-threaded warp operation.
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),
                             initializer=init_tile_worker,
                             initargs=(heights,validity,z,bounds,empty)) as pool:
        pending=set();iterator=iter(tasks);exhausted=False
        while pending or not exhausted:
            while not exhausted and len(pending)<2*workers:
                task=next(iterator,None)
                if task is None:exhausted=True
                else:pending.add(pool.submit(run_tile_task,task))
            if not pending:break
            done,pending=wait(pending,return_when=FIRST_COMPLETED)
            while done:
                future=done.pop()
                yield future.result()


def tiles(root,prefix='coastal',min_zoom=0,max_zoom=12,bounds=None,*,output_prefix=None,expected_inputs=None,workers=1):
    from build import atomic_json,digest_json,emit,read_json,sha256
    root=Path(root)
    if not isinstance(workers,int) or workers<1:raise ValueError('workers must be a positive integer')
    if not 0 <= min_zoom <= max_zoom <= 12:
        raise ValueError('Expected 0 <= min_zoom <= max_zoom <= 12')
    output_prefix=output_prefix or prefix
    height_path=root/f'mosaics/{prefix}.vrt'
    mask_path=root/f'mosaics/{prefix}-validity.vrt'
    recipe=dict(heights_vrt_sha256=sha256(height_path),validity_vrt_sha256=sha256(mask_path),
                code_sha256=sha256(Path(__file__)),encoder_sha256=sha256(Path(__file__).with_name('terrain_rgb.py')),
                gdal=rasterio.__gdal_version__,tile_size=TILE_SIZE,
                resampling='average z0-11, bilinear z12',validity='nearest pixel centre',
                warp_engine='GDAL chunked reproject per tile',warp_memory_limit_mib=32,
                low_zoom_crop='valid target bounding rectangle at z0-4; unchanged native source',
                minzoom=min_zoom,maxzoom=max_zoom,bounds=bounds,round_bits=QUANTIZATION,
                nodata_rgba=[1,134,160,255],nodata_preserved=False,river_flattening=False)
    river_path=root/f'rivers/{prefix}/cogs.json'
    if river_path.exists():
        river_manifest=read_json(river_path)
        river_source=river_manifest['river_source']
        if river_source not in ('glo30','wris'):
            raise ValueError('Unknown river mask in COG manifest')
        recipe.update(river_flattening=True,river_cog_manifest_sha256=sha256(river_path),
                      river_source=river_source,
                      river_scope='supported connected components; unsupported components unchanged',
                      quality_status=river_manifest.get('quality_status','Candidate river build; review before publication'))
    # The VRT bytes alone do not identify their dependent rasters.
    import xml.etree.ElementTree as ET
    sources=ET.parse(height_path).getroot().findall('.//SourceFilename')
    recipe['sheet_sha256']={Path(s.text).name:sha256(height_path.parent/s.text) for s in sources}
    masks=ET.parse(mask_path).getroot().findall('.//SourceFilename')
    recipe['validity_sha256']={Path(s.text).name:sha256(mask_path.parent/s.text) for s in masks}
    if expected_inputs is not None:
        # Extension must use the same rasters and decoding contract as the
        # completed zooms. Code/range/precision records are retained separately.
        for key in ('heights_vrt_sha256','validity_vrt_sha256','sheet_sha256','validity_sha256',
                    'encoder_sha256','gdal','tile_size','bounds','nodata_rgba','validity',
                    'river_flattening','river_cog_manifest_sha256','river_source'):
            if recipe.get(key)!=expected_inputs.get(key):
                raise ValueError(f'Cannot extend changed inputs or encoding: {key}')
    fingerprint=digest_json(recipe)
    tile_dir=root/'tiles'/output_prefix
    tile_dir.mkdir(parents=True,exist_ok=True)
    recipe_path=tile_dir/'recipe.json'
    if recipe_path.exists() and read_json(recipe_path)!=recipe:
        raise ValueError('Tile recipe or prepared inputs changed; use a new tile output prefix')
    atomic_json(recipe_path,recipe)
    empty=webp_bytes(np.tile(np.array([1,134,160,255],dtype=np.uint8),(TILE_SIZE,TILE_SIZE,1)))
    summaries=[]
    for z in range(min_zoom,max_zoom+1):
        started=time.monotonic()
        db=sqlite3.connect(tile_dir/f'z{z}.sqlite')
        db.execute('PRAGMA journal_mode=WAL')
        db.execute('PRAGMA synchronous=NORMAL')
        db.execute('CREATE TABLE IF NOT EXISTS tiles (tileid INTEGER PRIMARY KEY, x INTEGER, y INTEGER, data BLOB, sha256 TEXT, valid INTEGER, error REAL)')
        db.execute('CREATE TABLE IF NOT EXISTS batches (x INTEGER PRIMARY KEY, recipe TEXT)')
        db.execute('CREATE INDEX IF NOT EXISTS tiles_x ON tiles(x)')
        max_error=0.0
        completed=0
        with Sampler(height_path,mask_path,z,bounds=bounds) as sampler:
            g=sampler.grid
            expected_tiles=(g['xmax']-g['xmin']+1)*(g['ymax']-g['ymin']+1)
            pending_x=[]
            for x in range(g['xmin'],g['xmax']+1):
                done=db.execute('SELECT recipe FROM batches WHERE x=?',(x,)).fetchone()
                if done:
                    if done[0]!=fingerprint:
                        raise ValueError('Mismatching tile batch recipe')
                    row=db.execute('SELECT count(*),max(error) FROM tiles WHERE x=?',(x,)).fetchone()
                    if row[0]!=g['ymax']-g['ymin']+1:
                        raise ValueError('Incomplete committed tile batch')
                    completed+=row[0]
                    max_error=max(max_error,row[1] or 0)
                    continue
                pending_x.append(x)
            # A crash can leave completed tasks in an unfinished x-column.
            # Discard only those fragments; completed columns remain reusable.
            with db:
                for x in pending_x:db.execute('DELETE FROM tiles WHERE x=?',(x,))
            chunks_per_x=math.ceil((g['ymax']-g['ymin']+1)/ROWS_PER_TASK)
            remaining={x:chunks_per_x for x in pending_x}
            tasks=((x,y,min(g['ymax']+1,y+ROWS_PER_TASK)) for x in pending_x
                   for y in range(g['ymin'],g['ymax']+1,ROWS_PER_TASK))
            # Tiny global zooms scan most of the source; keep z0-4 serial to
            # avoid concurrent broad reads and their higher memory demand.
            used_workers=min(workers,len(pending_x)*chunks_per_x) if z>=5 else 1
            if used_workers>1:
                results=parallel_rows(height_path,mask_path,z,bounds,empty,tasks,used_workers)
            else:
                results=(encode_rows(sampler,empty,z,*task) for task in tasks)
            for x,rows in results:
                with db:
                    db.executemany('INSERT INTO tiles VALUES (?,?,?,?,?,?,?)',rows)
                    remaining[x]-=1
                    if remaining[x]==0:db.execute('INSERT INTO batches VALUES (?,?)',(x,fingerprint))
                completed+=len(rows)
                max_error=max(max_error,max(row[-1] for row in rows))
                state=dict(zoom=z,complete=completed,total=expected_tiles,
                           elapsed_seconds=round(time.monotonic()-started,1),max_encoding_error_m=max_error,
                           workers=used_workers)
                atomic_json(root/f'reports/{output_prefix}-tile-status.json',state)
                emit('tiles',**state)
        counts=db.execute('SELECT count(*),sum(length(data)),sum(valid),sum(valid=0) FROM tiles').fetchone()
        if counts[0]!=expected_tiles:
            raise ValueError('Tile coverage count differs')
        db.execute('PRAGMA wal_checkpoint(TRUNCATE)')
        db.close()
        report=dict(zoom=z,tiles=counts[0],bytes=counts[1],valid_pixels=counts[2],
                    fully_nodata_tiles=counts[3],max_encoding_error_m=max_error,
                    height_step_m=.1*2**round_bits(z),elapsed_seconds=round(time.monotonic()-started,2),
                    workers=used_workers)
        atomic_json(tile_dir/f'z{z}.json',report)
        summaries.append(report)
    atomic_json(root/f'reports/{output_prefix}-tiles.json',dict(recipe=recipe,zooms=summaries))
    return summaries


def package(root,prefix='coastal',max_bytes=1800000000):
    """Bounded-size PMTiles partitions, exact tile-byte verification, portable mosaic."""
    from build import atomic_json,emit,read_json,sha256
    from pmtiles.tile import Compression,TileType,tileid_to_zxy
    from pmtiles.reader import Reader,all_tiles
    from pmtiles.writer import Writer
    root=Path(root)
    tile_dir=root/'tiles'/prefix
    report=read_json(root/f'reports/{prefix}-tiles.json')
    recipe=report['recipe']
    outdir=root/'output'/prefix
    outdir.mkdir(parents=True,exist_ok=True)
    # PMTiles writer uses a temporary tile spool; keep it in the selected workspace.
    def local_writer(stream):
        previous=tempfile.tempdir
        try:
            tempfile.tempdir=str(root/'tmp')
            return Writer(stream)
        finally:
            tempfile.tempdir=previous
    metadata=dict(name='CartoDEM V3R1 EGM2008 TerrainRGB',format='webp',encoding='mapbox',
                  tileSize=TILE_SIZE,vertical_datum='EGM2008',alpha=255,nodata_height_m=0,nodata_preserved=False,
                  nodata_rgb=[1,134,160],quantization='round_bits='+recipe['round_bits'],river_flattening=recipe['river_flattening'],
                  coastal_mask='Copernicus GLO-30 ocean class 1',
                  attribution='CartoDEM: NRSC/ISRO; EGM2008: NGA; ocean mask: Copernicus DEM')
    if recipe['river_flattening']:
        description={'glo30':'Copernicus GLO-30 WBM class 3 (River); class 2 (Lake) excluded',
                     'wris':'WRIS River polygons; Sand and Island excluded'}[recipe['river_source']]
        metadata.update(river_mask=description,
                        quality_status=recipe['quality_status'])
        metadata['attribution']+=('; river mask: Copernicus GLO-30' if recipe['river_source']=='glo30'
                                  else '; river polygons: India-WRIS')
    partitions=[]
    writer=stream=None
    batch=[]
    byte_sum=0
    part=0
    def finish():
        nonlocal writer,stream,batch,byte_sum,part
        if not batch:
            return
        path=outdir/f'CartoDEM_v3r1_EGM2008_{part:03d}.pmtiles'
        temp=path.with_suffix('.pmtiles.part')
        zooms=[r[1] for r in batch]
        # Tile extent bounds describe actual archive coverage, including empty tiles.
        boxes=[]
        for tileid,z,x,y,digest in batch:
            n=2**z
            lon0=x/n*360-180;lon1=(x+1)/n*360-180
            lat1=math.degrees(math.atan(math.sinh(math.pi*(1-2*y/n))))
            lat0=math.degrees(math.atan(math.sinh(math.pi*(1-2*(y+1)/n))))
            boxes.append((lon0,lat0,lon1,lat1))
        bounds=[min(b[0] for b in boxes),min(b[1] for b in boxes),
                max(b[2] for b in boxes),max(b[3] for b in boxes)]
        header=dict(tile_type=TileType.WEBP,tile_compression=Compression.NONE,
                    min_lon_e7=round(bounds[0]*1e7),min_lat_e7=round(bounds[1]*1e7),
                    max_lon_e7=round(bounds[2]*1e7),max_lat_e7=round(bounds[3]*1e7),
                    min_zoom=min(zooms),max_zoom=max(zooms),center_zoom=min(zooms),
                    center_lon_e7=round((bounds[0]+bounds[2])*5e6),
                    center_lat_e7=round((bounds[1]+bounds[3])*5e6))
        meta={**metadata,'minzoom':min(zooms),'maxzoom':max(zooms),'bounds':bounds}
        writer.finalize(header,meta)
        stream.close()
        # Verify EVERY tile by checksum against the encoded tile database. Each
        # encoded image was also checked against its numeric elevations at creation.
        with temp.open('rb') as handle:
            def read_at(offset,length):
                handle.seek(offset)
                return handle.read(length)
            reader=Reader(read_at)
            actual=reader.header()
            if actual['tile_type']!=TileType.WEBP or actual['tile_compression']!=Compression.NONE:
                raise ValueError('Incorrect PMTiles format/compression header')
            expected={(z,x,y):digest for tileid,z,x,y,digest in batch}
            seen=set()
            for (z,x,y),data in all_tiles(read_at):
                if ((z,x,y) in seen or (z,x,y) not in expected
                        or hashlib.sha256(data).hexdigest()!=expected[(z,x,y)]):
                    raise ValueError(f'Packaged tile changed: {z}/{x}/{y}')
                seen.add((z,x,y))
            if len(seen)!=len(expected):
                raise ValueError('Packaged tile coverage differs')
        temp.replace(path)
        entry=dict(file=path.name,bytes=path.stat().st_size,sha256=sha256(path),
                   tiles=len(batch),verified_tiles=len(batch),bounds=bounds,
                   minzoom=min(zooms),maxzoom=max(zooms))
        partitions.append(entry)
        emit('package',**entry)
        writer=stream=None;batch=[];byte_sum=0;part+=1
    for zoom_report in report['zooms']:
        z=zoom_report['zoom']
        database=root/zoom_report.get('database',f'tiles/{prefix}/z{z}.sqlite')
        if not database.resolve().is_relative_to(root.resolve()):
            raise ValueError('Tile database must stay inside the build workspace')
        with closing(sqlite3.connect(f'file:{database}?mode=ro',uri=True)) as db:
            for tileid,x,y,data,digest in db.execute('SELECT tileid,x,y,data,sha256 FROM tiles ORDER BY tileid'):
                if hashlib.sha256(data).hexdigest()!=digest:
                    raise ValueError('Tile database checksum changed')
                if batch and byte_sum+len(data)>max_bytes:
                    finish()
                if writer is None:
                    temp=outdir/f'CartoDEM_v3r1_EGM2008_{part:03d}.pmtiles.part'
                    stream=temp.open('wb')
                    writer=local_writer(stream)
                writer.write_tile(tileid,data)
                batch.append((tileid,z,x,y,digest))
                byte_sum+=len(data)
    finish()
    atomic_json(outdir/'archives.json',dict(metadata=metadata,partitions=partitions,recipe=recipe,
                verified_tiles=sum(p['verified_tiles'] for p in partitions)))
    atomic_json(root/f'reports/{prefix}-package.json',dict(partitions=partitions,
                verified_tiles=sum(p['verified_tiles'] for p in partitions),
                total_bytes=sum(p['bytes'] for p in partitions)))
    return partitions
