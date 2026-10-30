"""Isolate per-sheet GDAL allocations while retaining global geometry IDs."""
import argparse
import math
from pathlib import Path
import time

import numpy as np
import rasterio
from rasterio.windows import Window

from build import atomic_json,atomic_npz,digest_json,read_json,sha256
from .geometry import extract_window
from .allocation import make_allocator


def geometry_sheet(root,work,sheet):
    started=time.monotonic();recipe=read_json(work/'geometry-recipe.json');halo=recipe['halo']
    with rasterio.Env(GDAL_CACHEMAX=48*1024**2,GDAL_NUM_THREADS='1',GDAL_MAX_DATASET_POOL_SIZE=24), \
         rasterio.open(work/'coastal.vrt') as heights,rasterio.open(work/'classes.vrt') as classes:
        t=rasterio.Affine(*sheet['transform'])
        x,y=(~heights.transform)*(t.c,t.f);x,y=round(x),round(y)
        x0=max(0,x-halo);y0=max(0,y-halo)
        x1=min(heights.width,x+sheet['width']+halo);y1=min(heights.height,y+sheet['height']+halo)
        window=Window(x0,y0,x1-x0,y1-y0)
        h=heights.read(1,window=window);mask=classes.read(1,window=window)
        core=(y-y0,y-y0+sheet['height'],x-x0,x-x0+sheet['width'])
        r0,r1,c0,c1=core
        if np.any((h[r0:r1,c0:c1]!=-255)&(mask[r0:r1,c0:c1]>3)):
            raise ValueError('Missing classification at valid elevation')
        lat=heights.transform.f+heights.transform.e*(y0+(y1-y0)/2)
        dy=abs(heights.transform.e)*111195.0802
        dx=heights.transform.a*111195.0802*math.cos(math.radians(lat))
        data=extract_window(h,mask,dx=dx,dy=dy,grid_width=heights.width,
                            row_offset=y0,col_offset=x0,core=core,
                            allocator=make_allocator(root,work,sheet['name'],recipe))
    path=work/'geometry'/(Path(sheet['name']).stem+'.npz');atomic_npz(path,**data)
    report=dict(sheet=sheet['name'],sha256=sha256(path),recipe_sha256=digest_json(recipe),
                river_pixels=len(data['pixels']),skeleton_pixels=len(data['nodes']),
                bank_samples=len(data['bank_heights']),elapsed_seconds=round(time.monotonic()-started,2))
    atomic_json(path.with_suffix('.json'),report)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--label',required=True)
    parser.add_argument('--sheet',required=True)
    args=parser.parse_args();root=args.workdir;work=root/'rivers'/args.label
    sheet=next(s for s in read_json(root/'reports/inputs.json')['sheets'] if s['name']==args.sheet)
    geometry_sheet(root,work,sheet)


if __name__=='__main__':main()
