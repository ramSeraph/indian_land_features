"""Evaluate ready-made Whitebox cost allocation without changing production caches."""
import argparse
import json
import math
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time

import numpy as np
import rasterio
from rasterio.windows import Window
from scipy.sparse.csgraph import dijkstra
from skimage.morphology import skeletonize

from build import atomic_json,read_json,sha256
from .geometry import pixel_graph,prune_spurs,extract_window
from .allocation import write_inputs as write_rasters,run_tools,read_owners,backend_recipe,make_allocator


def write_inputs(folder,river,seeds,dx,dy):
    folder.mkdir(parents=True,exist_ok=False)
    np.save(folder/'river.npy',river);np.save(folder/'seeds.npy',seeds)
    record=write_rasters(folder,river,seeds,dx,dy)
    record['files_sha256'].update({p:sha256(folder/p) for p in ('river.npy','seeds.npy')})
    atomic_json(folder/'inputs.json',record)


def run(root,folder):
    info=read_json(folder/'inputs.json')
    for p,digest in info['files_sha256'].items():
        if sha256(folder/p)!=digest:raise ValueError('Trial input changed')
    tool,results=run_tools(root,folder)
    river=np.load(folder/'river.npy',mmap_mode='r');seeds=np.load(folder/'seeds.npy',mmap_mode='r')
    owners=read_owners(folder,river,seeds)
    report=dict(state='complete',inputs=info,tool=tool,runs=results,assigned_river_pixels=len(owners),
                published_heights_changed=False,production_cache_changed=False,
                code_sha256=sha256(Path(__file__)),allocation=backend_recipe(root,'whitebox'),
                outputs_sha256={p:sha256(folder/p) for p in ('distance.tif','backlink.tif','owners.tif')})
    atomic_json(folder/'report.json',report)
    return report


def synthetic(root,folder):
    folder.mkdir(parents=True,exist_ok=False);results=[]
    for name,dx,dy in [('island',23.,31.),('parallel',23.,31.),('ties',30.,30.),('diagonal',23.,31.)]:
        river=np.zeros((81,111),bool);seeds=np.zeros_like(river)
        if name=='diagonal':
            river[np.arange(5,76),np.arange(15,86)]=True
            seeds[5,15]=True;seeds[75,85]=True
        elif name=='parallel':
            river[8:20,5:100]=True;river[40:55,5:100]=True
            seeds[14,10]=True;seeds[47,90]=True
        else:
            river[5:76,5:106]=True
            if name=='island':river[20:62,30:82]=False
            seeds[10,10]=True;seeds[70,100]=True
        case=folder/name;write_inputs(case,river,seeds,dx,dy);report=run(root,case)
        flat=np.flatnonzero(river);seed_pixels=np.flatnonzero(seeds)
        graph=pixel_graph(flat,river.shape[1],dx,dy)
        indices=np.searchsorted(flat,seed_pixels)
        expected=dijkstra(graph,directed=False,indices=indices)
        nearest=expected.min(axis=0)
        with rasterio.open(case/'distance.tif') as ds:actual=ds.read(1).ravel()[flat]
        with rasterio.open(case/'owners.tif') as ds:labels=ds.read(1).ravel()[flat].astype(int)-1
        np.testing.assert_allclose(actual,nearest,rtol=2e-7,atol=1e-5)
        np.testing.assert_allclose(expected[labels,np.arange(len(flat))],nearest,rtol=1e-12,atol=1e-8)
        _,_,old=dijkstra(graph,directed=False,indices=indices,min_only=True,return_predecessors=True)
        different=int(np.sum(seed_pixels[labels]!=flat[old]))
        results.append(dict(case=name,river_pixels=len(flat),maximum_distance_rounding_m=float(np.max(np.abs(actual-nearest))),
                            owner_differences=different,all_assigned_sources_minimum_cost=True,
                            mask_preserved=True,runs=report['runs']))
    # Exercise the exact production adapter, including its scratch cleanup and
    # shared global pixel IDs across two overlapping processing windows.
    recipe=dict(allocation=backend_recipe(root,'whitebox'))
    classes=np.zeros((50,80),dtype=np.uint8);classes[21:28,3:77]=3
    heights=np.full(classes.shape,10.,dtype=np.float32)
    pieces=[]
    for name,x0,x1,core in [('whole',0,80,(0,50,0,80)),('west',0,56,(0,50,0,40)),
                           ('east',24,80,(0,50,16,56))]:
        pieces.append(extract_window(heights[:,x0:x1],classes[:,x0:x1],dx=23,dy=31,
            grid_width=80,col_offset=x0,core=core,
            allocator=make_allocator(root,folder,name,recipe)))
    whole,west,east=pieces
    pixels=np.r_[west['pixels'],east['pixels']];order=np.argsort(pixels)
    np.testing.assert_array_equal(pixels[order],whole['pixels'])
    np.testing.assert_array_equal(np.r_[west['owners'],east['owners']][order],whole['owners'])
    np.testing.assert_array_equal(np.sort(np.r_[west['nodes'],east['nodes']]),whole['nodes'])
    result=dict(cases=results,production_adapter_cross_window_check=True)
    atomic_json(folder/'summary.json',result);print(json.dumps(result),flush=True)


def prepare(root,folder,sheet,label):
    work=root/'rivers'/label;grid=read_json(work/'grid.json')
    item=next(r for r in read_json(root/'reports/inputs.json')['sheets'] if r['name']==sheet)
    halo=read_json(work/'geometry-recipe.json')['halo'];native=rasterio.Affine(*grid['transform'])
    t=rasterio.Affine(*item['transform']);x,y=(~native)*(t.c,t.f);x,y=round(x),round(y)
    x0,y0=max(0,x-halo),max(0,y-halo)
    x1,y1=min(grid['width'],x+item['width']+halo),min(grid['height'],y+item['height']+halo)
    with rasterio.Env(GDAL_CACHEMAX=32*1024**2),rasterio.open(work/'classes.vrt') as ds:
        river=ds.read(1,window=Window(x0,y0,x1-x0,y1-y0))==3
    lat=native.f+native.e*(y0+(y1-y0)/2)
    dy=abs(native.e)*111195.0802;dx=native.a*111195.0802*math.cos(math.radians(lat))
    seeds,pruned=prune_spurs(skeletonize(river),river,dx,dy)
    write_inputs(folder,river,seeds,dx,dy)
    atomic_json(folder/'context.json',dict(sheet=sheet,label=label,window=[x0,y0,x1-x0,y1-y0],
                pruned_skeleton_pixels=pruned,geometry_code_sha256=sha256(Path(__file__).with_name('geometry.py')),
                peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
    print(json.dumps(read_json(folder/'inputs.json')),flush=True)


def worker(root,folder,sheet,label):
    """Benchmark the real isolated worker on an unchanged existing mosaic context."""
    previous=root/'rivers'/label
    folder.mkdir(parents=True,exist_ok=False)
    (folder/'geometry').mkdir()
    for name in ('coastal.vrt','classes.vrt','validity.vrt','grid.json'):
        shutil.copy2(previous/name,folder/name)
    recipe=read_json(previous/'geometry-recipe.json')
    recipe.update(allocation=backend_recipe(root,'whitebox'),
                  benchmark_context=dict(label=label,recipe_sha256=sha256(previous/'geometry-recipe.json')),
                  geometry_sha256=sha256(Path(__file__).with_name('geometry.py')),
                  staging_sha256=sha256(Path(__file__).with_name('staging.py')))
    atomic_json(folder/'geometry-recipe.json',recipe)
    started=time.monotonic();peak=0;individual=0;stopped=False
    with (folder/'worker.log').open('w') as log:
        process=subprocess.Popen([sys.executable,'-m','river_flattening.staging','geometry',
            '--workdir',str(root),'--label',folder.name,'--sheet',sheet],
            stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        try:
            while process.poll() is None:
                rows=[list(map(int,line.split())) for line in subprocess.check_output(
                    ['ps','-eo','pid=,ppid=,rss='],text=True).splitlines()]
                family={process.pid}
                while True:
                    children={pid for pid,parent,rss in rows if parent in family}
                    if children<=family:break
                    family.update(children)
                rss=[rss for pid,parent,rss in rows if pid in family]
                peak=max(peak,sum(rss));individual=max(individual,max(rss,default=0))
                if sum(rss)>1700*1024:
                    stopped=True;os.killpg(process.pid,signal.SIGTERM);break
                time.sleep(.05)
            code=process.wait(timeout=10)
        finally:
            if process.poll() is None:
                os.killpg(process.pid,signal.SIGKILL);process.wait()
    report=dict(sheet=sheet,source_context=label,exit_code=code,
        elapsed_seconds=round(time.monotonic()-started,2),peak_observed_tree_rss_mib=round(peak/1024,2),
        peak_observed_individual_rss_mib=round(individual/1024,2),sampling_interval_seconds=.05,
        stopped_for_memory=stopped,production_cache_changed=False,
        code_sha256=sha256(Path(__file__)))
    atomic_json(folder/'worker-report.json',report);print(json.dumps(report),flush=True)
    if code:raise RuntimeError('Production geometry worker trial failed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['synthetic','prepare','run','worker'])
    parser.add_argument('--workdir',required=True,type=Path)
    parser.add_argument('--label',default='allocation-trial-v1')
    parser.add_argument('--geometry-label',default='final')
    parser.add_argument('--sheet',default='cdnf46g.tif')
    args=parser.parse_args();root=args.workdir.resolve()
    if Path(args.label).name!=args.label:parser.error('Use a plain label')
    folder=root/'rivers'/args.label
    if args.stage=='synthetic':synthetic(root,folder)
    elif args.stage=='prepare':prepare(root,folder,args.sheet,args.geometry_label)
    elif args.stage=='worker':worker(root,folder,args.sheet,args.geometry_label)
    else:print(json.dumps(run(root,folder)),flush=True)


if __name__=='__main__':main()
