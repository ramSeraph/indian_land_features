"""GRASS routing evidence on a disposable surface, never on published heights.

This bounded regional runner uses the native geographic grid. It interpolates
small internal holes only for routing. Sea stays NULL; original validity is
retained in a separate evidence raster. Output elevation COGs never consume the
interpolated raster. D8 drainage is supporting evidence, not a river mask.
"""
import os
from pathlib import Path
import shutil
import subprocess
import sys

if __package__ in (None, ''):
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import numpy as np
import rasterio
from rasterio.fill import fillnodata
from rasterio.windows import Window
from scipy.ndimage import binary_dilation

from build import atomic_json, read_json, sha256

# GRASS positive directions, counterclockwise from northeast.
OFFSETS = np.array([[0,0],[-1,1],[-1,0],[-1,-1],[0,-1],
                    [1,-1],[1,0],[1,1],[0,1]],dtype=np.int8)


def analysis_window(height, classes, max_search=64):
    valid = np.isfinite(height)&(height!=-255)&(classes!=1)&(classes<4)
    h = np.where(valid,height,0).astype(np.float32)
    if valid.any():
        h = fillnodata(h, mask=valid.astype('uint8'), max_search_distance=max_search)
        # GDAL leaves distant holes at zero: identify fillable cells separately.
        from scipy.ndimage import distance_transform_edt
        filled = (~valid)&(classes!=1)&(classes<4)&(distance_transform_edt(~valid)<=max_search)
    else:
        filled = np.zeros(valid.shape,bool)
    h[~(valid|filled)] = -255
    # Exclude water-flow samples near either provisional interpolation or holes.
    uncertain = binary_dilation((~valid)&(classes!=1), iterations=8)
    evidence = valid&~uncertain
    return h,evidence,filled


def grass_worker(directory):
    """Invoked by GRASS --exec using its isolated mapset."""
    d=Path(directory).resolve()
    commands=[
        ['r.in.gdal',f'input={d/"analysis.tif"}','output=elevation','--overwrite','--quiet'],
        ['g.region','raster=elevation'],
        ['r.watershed','-sm','elevation=elevation','drainage=drainage','accumulation=accumulation',
         'memory=256','--overwrite','--quiet'],
    ]
    for name,dtype,nodata in [('drainage','Int16',-32768),('accumulation','Float64',-999999)]:
        commands.append(['r.out.gdal',f'input={name}',f'output={d/(name+".tif")}',
                         f'type={dtype}',f'nodata={nodata}','createopt=COMPRESS=DEFLATE,TILED=YES',
                         '--overwrite','--quiet'])
    for command in commands:
        subprocess.run(command,check=True)


def route(root,label,max_cells=40_000_000):
    work=Path(root)/'rivers'/label
    grid=read_json(work/'grid.json')
    if grid['width']*grid['height']>max_cells:
        from .regional import route_regions
        return route_regions(root,label,max_cells=max_cells)
    output=work/'routing'
    output.mkdir(exist_ok=True)
    grass=shutil.which('grass')
    if not grass:
        raise ValueError('Install the GRASS CLI before the routing stage')
    version=subprocess.check_output([grass,'--config','version'],text=True).strip()
    recipe=dict(geometry=read_json(work/'geometry.json')['recipe'],grass_version=version,
                code_sha256=sha256(Path(__file__)),algorithm='r.watershed -sm; memory=256',
                analysis_fill_distance_pixels=64, evidence_exclusion_pixels=8)
    report_path=output/'report.json'
    if report_path.exists():
        old=read_json(report_path)
        if old['recipe']!=recipe or any(sha256(output/p)!=h for p,h in old['outputs'].items()):
            raise ValueError('Routing recipe or output changed; use a new label')
        return old
    filled_total=0
    with rasterio.Env(GDAL_CACHEMAX=64*1024**2):
        with rasterio.open(work/'coastal.vrt') as src,rasterio.open(work/'classes.vrt') as mask:
            if src.width*src.height>max_cells:
                raise ValueError('Route bounded regions with overlapping context; this pilot runner is limited to 40 million cells')
            profile=dict(driver='GTiff',width=src.width,height=src.height,transform=src.transform,crs=src.crs,
                         count=1,tiled=True,blockxsize=256,blockysize=256,compress='DEFLATE')
            with rasterio.open(output/'analysis.tif','w',**profile,dtype='float32',nodata=-255,predictor=3) as dst, \
                 rasterio.open(output/'evidence.tif','w',**profile,dtype='uint8') as ev:
                for y in range(0,src.height,512):
                    for x in range(0,src.width,512):
                        w=min(512,src.width-x);h=min(512,src.height-y)
                        x0=max(0,x-72);y0=max(0,y-72)
                        x1=min(src.width,x+w+72);y1=min(src.height,y+h+72)
                        window=Window(x0,y0,x1-x0,y1-y0)
                        data,okay,filled=analysis_window(src.read(1,window=window),mask.read(1,window=window))
                        crop=np.s_[y-y0:y-y0+h,x-x0:x-x0+w]
                        dst.write(data[crop],1,window=Window(x,y,w,h))
                        ev.write(okay[crop].astype('uint8'),1,window=Window(x,y,w,h))
                        filled_total+=int(filled[crop].sum())
    # No shell interpolation, and no absolute execution paths in the report.
    env=os.environ.copy();env['GRASS_MESSAGE_FORMAT']='plain';env['OMP_NUM_THREADS']='1'
    with (output/'grass.log').open('w') as log:
        subprocess.run([grass,'--tmp-location',str((output/'analysis.tif').resolve()),'--exec',
                        sys.executable,str(Path(__file__).resolve()),'--worker',str(output.resolve())],
                       env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    report=dict(recipe=recipe,analysis_filled_pixels=filled_total,
                outputs={p:sha256(output/p) for p in ('analysis.tif','evidence.tif','drainage.tif','accumulation.tif')},
                published_elevations_interpolated=False)
    atomic_json(report_path,report)
    return report


def sample_owners(work,native,width):
    """Accumulate native flow-vector evidence at in-water skeleton owners.

    Positive accumulation only: GRASS flags incomplete boundary catchments with
    a negative value. Missing/interpolated heights are never flow observations.
    Weights are log(1+area), capped to avoid one channel pixel dominating a reach.
    """
    from .production import atomic_npz
    work=Path(work);directory=work/'routing'
    nodes=native['nodes']; vector=np.zeros((len(nodes),2));weight=np.zeros(len(nodes));counts=np.zeros(len(nodes),int)
    with rasterio.open(directory/'drainage.tif') as ds,rasterio.open(directory/'accumulation.tif') as ac, \
         rasterio.open(directory/'evidence.tif') as ev:
        records=read_json(work/'geometry.json')['sheets']
        for record in records:
            with np.load(work/'geometry'/(Path(record['sheet']).stem+'.npz')) as data:
                pixels,owners=data['pixels'],data['owners']
            if not len(pixels):
                continue
            rows,cols=pixels//width,pixels%width
            y,x=int(rows.min()),int(cols.min());h=int(rows.max()-y+1);w=int(cols.max()-x+1)
            window=Window(x,y,w,h);r,c=rows-y,cols-x
            direction=ds.read(1,window=window)[r,c];area=ac.read(1,window=window)[r,c]
            okay=(ev.read(1,window=window)[r,c]>0)&(direction>0)&(direction<=8)&(area>0)
            index=np.searchsorted(nodes,owners[okay])
            if np.any(index>=len(nodes)) or np.any(nodes[index]!=owners[okay]):
                raise ValueError('A routing sample has no shared skeleton owner')
            wt=np.minimum(np.log1p(area[okay]),12.)
            offsets=OFFSETS[direction[okay]]
            for axis in range(2):
                np.add.at(vector[:,axis],index,wt*offsets[:,axis])
            np.add.at(weight,index,wt);np.add.at(counts,index,1)
    atomic_npz(directory/'owner-evidence.npz',nodes=nodes,vector=vector,weight=weight,count=counts)
    return dict(vector=vector,weight=weight,count=counts)


if __name__=='__main__' and len(sys.argv)==3 and sys.argv[1]=='--worker':
    grass_worker(sys.argv[2])
