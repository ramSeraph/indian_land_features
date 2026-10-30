"""Bounded overlapping routing contexts; only reliable direction samples survive.

Each sheet owns its samples once. Artificial context edges never become claimed
river outlets. Incomplete upstream terrain is excluded from the evidence.
The globally connected river fit remains independent of these processing windows.
"""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
import rasterio
from rasterio.windows import Window

from build import atomic_json, digest_json, emit, read_json, sha256
from .production import atomic_npz
from .routing import analysis_window, OFFSETS


def accumulate(owners, direction, area, valid):
    good=valid&(direction>0)&(direction<=8)&(area>0)&np.isfinite(area)
    ids,index=np.unique(owners[good],return_inverse=True)
    vector=np.zeros((len(ids),2));weight=np.zeros(len(ids));count=np.zeros(len(ids),np.int64)
    wt=np.minimum(np.log1p(area[good]),12.)
    for axis in range(2):
        np.add.at(vector[:,axis],index,wt*OFFSETS[direction[good],axis])
    np.add.at(weight,index,wt);np.add.at(count,index,1)
    return dict(nodes=ids,vector=vector,weight=weight,count=count)


def route_regions(root,label,halo=1024,max_cells=40_000_000,backend='whitebox'):
    root=Path(root);work=root/'rivers'/label;output=work/'routing'
    output.mkdir(exist_ok=True)
    grass=None
    if backend=='whitebox':
        from .whitebox import installation,route_context
        tool,_,_=installation(root)
        algorithm='Whitebox BreachDepressions + D8Pointer + D8FlowAccumulation --pntr; cells'
    elif backend=='grass':
        grass=shutil.which('grass')
        if not grass:raise ValueError('Install GRASS for the historical GRASS backend')
        tool=dict(version=subprocess.check_output([grass,'--config','version'],text=True).strip())
        algorithm='r.watershed -sm; memory=256'
    else:raise ValueError('Unknown routing backend')
    records=read_json(work/'geometry.json');grid=read_json(work/'grid.json')
    transform=rasterio.Affine(*grid['transform'])
    recipe=dict(geometry_sha256=sha256(work/'geometry.json'),
                code={p:sha256(Path(__file__).with_name(p)) for p in ('regional.py','routing.py','whitebox.py')},
                backend=backend,tool=tool,
                halo_pixels=halo,max_cells=max_cells,algorithm=algorithm,
                only_positive_complete_accumulations=True)
    report_path=output/'report.json';recipe_id=digest_json(recipe)
    if report_path.exists():
        old=read_json(report_path)
        if old['recipe']!=recipe:raise ValueError('Regional routing recipe changed; use a new label')
        for r in old['regions']:
            if sha256(output/r['evidence'])!=r['sha256']:raise ValueError('Regional evidence changed')
        return old
    reports=[]
    env=os.environ.copy();env.update(GRASS_MESSAGE_FORMAT='plain',OMP_NUM_THREADS='1',TMPDIR=str((root/'tmp').resolve()))
    with rasterio.Env(GDAL_CACHEMAX=48*1024**2,GDAL_NUM_THREADS='1',GDAL_MAX_DATASET_POOL_SIZE=24), \
         rasterio.open(work/'coastal.vrt') as src,rasterio.open(work/'classes.vrt') as masks:
        for record in records['sheets']:
            name=record['sheet'];stem=Path(name).stem
            target=output/(stem+'-owners.npz');done=target.with_suffix('.json')
            if target.exists():
                report=read_json(done)
                if report['recipe_sha256']!=recipe_id or report['sha256']!=sha256(target):
                    raise ValueError('Regional routing cache differs')
                reports.append(report);continue
            started=time.monotonic()
            with np.load(work/'geometry'/(stem+'.npz')) as data:
                pixels,owners=data['pixels'],data['owners']
            filled_total=0;rasters={};bounds=None;routing_details=None
            if not len(pixels):
                evidence=accumulate(owners,np.empty(0,np.int16),np.empty(0),np.empty(0,bool))
            else:
                with rasterio.open(root/'sheets/coastal'/name) as sheet:
                    x,y=(~transform)*(sheet.transform.c,sheet.transform.f);x,y=round(x),round(y)
                    x0=max(0,x-halo);y0=max(0,y-halo)
                    x1=min(src.width,x+sheet.width+halo);y1=min(src.height,y+sheet.height+halo)
                width,height=x1-x0,y1-y0
                if width*height>max_cells:raise ValueError('Regional context exceeds configured cell budget')
                bounds=[x0,y0,width,height]
                temporary=root/'tmp'/('routing-'+label+'-'+stem)
                if temporary.exists():shutil.rmtree(temporary)
                temporary.mkdir()
                profile=dict(driver='GTiff',width=width,height=height,crs=src.crs,
                             transform=src.window_transform(Window(*bounds)),count=1,
                             tiled=True,blockxsize=256,blockysize=256,compress='DEFLATE')
                completed=False
                try:
                    with rasterio.open(temporary/'analysis.tif','w',**profile,dtype='float32',nodata=-255,predictor=3) as dst, \
                         rasterio.open(temporary/'evidence.tif','w',**profile,dtype='uint8') as ev:
                        for yy in range(0,height,512):
                            for xx in range(0,width,512):
                                w=min(512,width-xx);h=min(512,height-yy)
                                a=max(0,xx-72);b=max(0,yy-72)
                                c=min(width,xx+w+72);d=min(height,yy+h+72)
                                window=Window(x0+a,y0+b,c-a,d-b)
                                values,valid,filled=analysis_window(src.read(1,window=window),masks.read(1,window=window))
                                crop=np.s_[yy-b:yy-b+h,xx-a:xx-a+w]
                                dst.write(values[crop],1,window=Window(xx,yy,w,h))
                                ev.write(valid[crop].astype('uint8'),1,window=Window(xx,yy,w,h))
                                filled_total+=int(filled[crop].sum())
                    if backend=='whitebox':
                        routing_details=route_context(root,temporary,output,stem+'-whitebox-')
                        validity_name='routing-valid'
                    else:
                        with (output/(stem+'-grass.log')).open('w') as log:
                            subprocess.run([grass,'--tmp-location',str((temporary/'analysis.tif').resolve()),'--exec',
                                            sys.executable,str(Path(__file__).with_name('routing.py').resolve()),
                                            '--worker',str(temporary.resolve())],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
                        validity_name='evidence'
                    rows,cols=pixels//grid['width']-y0,pixels%grid['width']-x0
                    yy,xx=int(rows.min()),int(cols.min())
                    window=Window(xx,yy,int(cols.max()-xx+1),int(rows.max()-yy+1))
                    arrays=[]
                    for basename in ('drainage','accumulation',validity_name):
                        path=temporary/(basename+'.tif');rasters[path.name]=sha256(path)
                        with rasterio.open(path) as ds:arrays.append(ds.read(1,window=window)[rows-yy,cols-xx])
                    direction,area,valid=arrays
                    # Never use vectors adjacent to the artificial context edge.
                    valid=(valid>0)&(rows>=8)&(cols>=8)&(rows<height-8)&(cols<width-8)
                    evidence=accumulate(owners,direction,area,valid)
                    completed=True
                finally:
                    if completed:shutil.rmtree(temporary)
            atomic_npz(target,**evidence)
            report=dict(sheet=name,evidence=target.name,sha256=sha256(target),recipe_sha256=recipe_id,
                        context_window=bounds,analysis_filled_pixels=filled_total,
                        sampled_pixels=int(evidence['count'].sum()),owners=len(evidence['nodes']),
                        routing_details=routing_details,
                        temporary_raster_sha256=rasters,elapsed_seconds=round(time.monotonic()-started,2))
            atomic_json(done,report);reports.append(report)
            emit('regional_routing',complete=len(reports),total=len(records['sheets']),**report)
            atomic_json(output/'status.json',dict(complete=len(reports),total=len(records['sheets']),sheet=name))
    report=dict(mode='regional',recipe=recipe,regions=reports,published_elevations_interpolated=False)
    atomic_json(report_path,report)
    return report
