#!/usr/bin/env python3
"""Rebuild opaque TerrainRGB from existing final COGs, without DEM processing.

Use a fresh --workdir for the display rebuild. --cog-manifest is the original
rivers/final/cogs.json; every COG is checked against it before being used.
"""
import argparse
from contextlib import ExitStack
from datetime import datetime,timezone
import fcntl
import os
from pathlib import Path
import resource
import subprocess
import sys
import tarfile
import traceback

from build import atomic_json,create_mosaic,emit,read_json,setup,sha256
from validation.check_tiles_environment import check
from production_tiles import tiles,package


def rebuild(root,cog_dir,cog_manifest,*,label='display',expected_sheets=543,min_zoom=0,max_zoom=12,workers=1,tiles_only=True):
    root=setup(Path(root).resolve());cog_dir=Path(cog_dir).resolve()
    source_manifest=read_json(cog_manifest)
    entries=source_manifest['sheets']
    if len(entries)!=expected_sheets:
        raise ValueError(f'Expected {expected_sheets} COGs, found {len(entries)}')
    if source_manifest.get('river_source') not in ('glo30','wris'):
        raise ValueError('Expected the final river-treated COG manifest')
    names=[e['file'] for e in entries]
    if len(set(names))!=len(names) or any(Path(n).name!=n or not n.endswith('.tif') for n in names):
        raise ValueError('COG filenames must be unique plain TIFF names')
    paths=[];inputs=[]
    for entry in sorted(entries,key=lambda e:e['file']):
        path=cog_dir/entry['file'];digest=sha256(path)
        if digest!=entry['sha256'] or path.stat().st_size!=entry['bytes']:
            raise ValueError(f'COG differs from its original manifest: {path.name}')
        paths.append(path)
        inputs.append(dict(file=path.name,bytes=entry['bytes'],sha256=digest))
    record=dict(cog_directory=os.path.relpath(cog_dir,root),
                original_manifest_sha256=sha256(cog_manifest),sheets=inputs)
    input_path=root/f'reports/{label}-cog-inputs.json'
    if input_path.exists() and read_json(input_path)!=record:
        raise ValueError('COG inputs changed; use a fresh output workspace')
    atomic_json(input_path,record)
    atomic_json(root/f'rivers/{label}/cogs.json',source_manifest)
    emit('verified_cogs',sheets=len(inputs),bytes=sum(e['bytes'] for e in inputs))

    height=root/f'mosaics/{label}.vrt';mask=root/f'mosaics/{label}-validity.vrt'
    mosaic_path=root/f'reports/{label}-mosaic.json'
    if mosaic_path.exists():
        previous=read_json(mosaic_path)
        for path,key in ((height,'height_sha256'),(mask,'mask_sha256')):
            if sha256(path)!=previous[key]:raise ValueError('Prepared mosaic changed')
    else:
        emit('mosaic_start')
        report=create_mosaic(paths,height,mask)
        report.update(height_sha256=sha256(height),mask_sha256=sha256(mask))
        atomic_json(mosaic_path,report)
    emit('tiles_start',minzoom=min_zoom,maxzoom=max_zoom)
    tile_stats=tiles(root,label,min_zoom,max_zoom,workers=workers)
    if tiles_only:
        # The final packaging workflow exports these databases to MBTiles and
        # invokes the standard pmtiles-mosaic CLI. Do not run the legacy writer.
        partitions=[]
        recipe=read_json(root/f'reports/{label}-tiles.json')['recipe']
        metadata=dict(name='CartoDEM V3R1 EGM2008 TerrainRGB',format='webp',encoding='mapbox',
                      tileSize=512,vertical_datum='EGM2008',alpha=255,nodata_height_m=0,
                      nodata_preserved=False,nodata_rgb=[1,134,160],quantization='round_bits=16-max(z,5)',
                      attribution='CartoDEM: NRSC/ISRO; EGM2008: NGA; ocean/river masks: Copernicus GLO-30',
                      quality_status=source_manifest.get('quality_status',''))
        metadata.update(river_flattening=recipe['river_flattening'],
                        coastal_mask='Copernicus GLO-30 ocean class 1',
                        river_mask={'glo30':'Copernicus GLO-30 WBM class 3 (River); class 2 (Lake) excluded',
                                    'wris':'WRIS River polygons; Sand and Island excluded'}[recipe['river_source']])
        atomic_json(root/f'output/{label}/archives.json',dict(metadata=metadata,partitions=[],
                                                           recipe=recipe,verified_tiles=0))
    else:
        emit('package_start')
        partitions=package(root,label,max_bytes=1_800_000_000)
    # Independent checksums for the COG release; files remain on their original disk.
    release_dir=root/f'output/{label}'
    atomic_json(release_dir/'cogs.json',source_manifest)
    (release_dir/'COG_SHA256SUMS.txt').write_text(''.join(f"{e['sha256']}  {e['file']}\n" for e in inputs))
    (release_dir/'PMTILES_SHA256SUMS.txt').write_text(''.join(f"{e['sha256']}  {e['file']}\n" for e in partitions))
    software=check();atomic_json(release_dir/'software.json',software)
    code=Path(__file__).resolve().parent
    source_names=('rebuild_tiles.py','build.py','cartodem_to_egm2008.py','production_tiles.py',
                  'terrain_rgb.py','validation/__init__.py','validation/check_tiles_environment.py','requirements-tiles.txt',
                  'install_tiles_arm64.sh','TILING_SETUP.md',
                  'export_mbtiles.py','partition_mbtiles.sh','validation/verify_partitioned.py',
                  'validation/verify_mosaic.mjs','requirements-mosaic.txt')
    with tarfile.open(release_dir/'tiling-source.tar.gz','w:gz') as archive:
        for name in source_names:archive.add(code/name,arcname=name,recursive=False)
    notes=f'''# CartoDEM display tiles rebuilt from final COGs

Completed {datetime.now(timezone.utc).isoformat()}.

The {len(inputs)} input COGs were checked against their original manifest
({record['original_manifest_sha256']}). The COGs are unchanged: native Float32
EGM2008 elevations and declared NoData=-255, after ocean and river treatment.
The original treatment metadata and quality status are retained in cogs.json.
No datum conversion, water treatment or routing was repeated.

Display tiles use Mapbox TerrainRGB, 512 pixels, zooms {min_zoom}–{max_zoom},
lossless WebP and round_bits=16-max(z,5). The z12 increment is 1.6 m, with
at most 0.8 m rounding error. Missing pixels become RGB (1,134,160), exactly
0 metres. Every pixel is opaque (alpha 255); display tiles carry no nodata
marker. Numeric resampling excludes missing samples before the placeholder
is inserted. Native validity is resampled separately by nearest pixel centre.

Tile databases resume at completed x-column batches. Each encoded image is
decoded and checked byte for byte, and every packaged tile is checked against
its database. {'Export the databases to MBTiles and run pmtiles-mosaic partition with a 2,000,000,000-byte limit.' if tiles_only else 'This historical packaging path targets 1.8 billion payload bytes; use --tiles-only and the standard pmtiles-mosaic workflow for publication.'}
The requested worker count was {workers}. Zooms 0–4 and packaging run serially;
other zooms use isolated processes with one GDAL thread each. Worker count does
not alter the encoding recipe and can be changed when resuming.
COG and PMTiles checksums are recorded separately. software.json records the
runtime; tiling-source.tar.gz contains the rebuild code and instructions.

| Archive | Bytes | SHA-256 |
| --- | --- | --- |
'''
    for entry in partitions:
        notes+=f"| {entry['file']} | {entry['bytes']} | `{entry['sha256']}` |\n"
    (release_dir/'TILE_BUILD_NOTES.md').write_text(notes)
    return dict(sheets=len(inputs),tiles=sum(e['tiles'] for e in tile_stats),
                archives=len(partitions),total_bytes=sum(e['bytes'] for e in partitions))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',required=True,type=Path)
    parser.add_argument('--cog-dir',required=True,type=Path)
    parser.add_argument('--cog-manifest',required=True,type=Path)
    parser.add_argument('--label',default='display')
    parser.add_argument('--expected-sheets',type=int,default=543)
    parser.add_argument('--min-zoom',type=int,default=0)
    parser.add_argument('--max-zoom',type=int,default=12)
    parser.add_argument('--workers',type=int,default=1,help='Parallel tile processes at z5-12 (default: 1)')
    parser.add_argument('--tiles-only',action='store_true',default=True,
                        help='Leave encoded databases for MBTiles export and pmtiles-mosaic (default)')
    parser.add_argument('--detach',action='store_true')
    args=parser.parse_args()
    if not args.label or Path(args.label).name!=args.label or args.label in ('.','..'):
        parser.error('Use a plain output label')
    if not 0<=args.min_zoom<=args.max_zoom<=12:parser.error('Zooms must be between 0 and 12')
    if args.expected_sheets<1:parser.error('Expected sheet count must be positive')
    if args.workers<1:parser.error('Workers must be positive')
    root=setup(args.workdir.expanduser().resolve())
    # Refuse the existing output tree: legacy tiles must not be mixed with opaque tiles.
    recipe=root/f'tiles/{args.label}/recipe.json'
    if recipe.exists() and read_json(recipe).get('nodata_rgba')!=[1,134,160,255]:
        parser.error('Existing tile output uses a different nodata encoding; select a fresh workdir')
    if args.detach:
        command=[sys.executable,str(Path(__file__).resolve()),*[a for a in sys.argv[1:] if a!='--detach']]
        with (root/f'logs/rebuild-{args.label}.log').open('a') as log:
            child=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=log,stderr=log,start_new_session=True)
        emit('started',pid=child.pid,status=f'reports/rebuild-{args.label}.json');return
    record=dict(state='running',pid=os.getpid(),started=datetime.now(timezone.utc).isoformat(),
                code_sha256=sha256(Path(__file__)),workers=args.workers,tiles_only=args.tiles_only)
    status=root/f'reports/rebuild-{args.label}.json'
    with ExitStack() as stack:
        for name in ('build.lock',f'run-{args.label}.lock'):
            lock=stack.enter_context((root/'logs'/name).open('a'))
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        atomic_json(status,record)
        try:
            result=rebuild(root,args.cog_dir,args.cog_manifest,label=args.label,
                           expected_sheets=args.expected_sheets,min_zoom=args.min_zoom,max_zoom=args.max_zoom,
                           workers=args.workers,tiles_only=args.tiles_only)
            record.update(state='complete',**result)
        except BaseException as error:
            record.update(state='failed',error=str(error));traceback.print_exc();raise
        finally:
            record['finished']=datetime.now(timezone.utc).isoformat()
            record['peak_rss_kib']=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024 if sys.platform=='darwin' else 1))
            record['rss_scope']='parent process; excludes concurrent workers'
            record['child_ru_maxrss_kib']=int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/(1024 if sys.platform=='darwin' else 1))
            atomic_json(status,record)


if __name__=='__main__':main()
