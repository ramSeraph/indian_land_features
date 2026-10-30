"""Prepare GLO-30 river geometry and bank observations on the shared native grid."""
import argparse
import importlib.metadata
import inspect
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np
import rasterio
from rasterio.windows import Window

from build import atomic_json,atomic_npz,create_mosaic,digest_json,emit,read_json,sha256,setup
from .geometry import extract_window
from .allocation import allocation_recipe,make_allocator


def mask_mosaic(root,selected,report,path):
    transform=rasterio.Affine(*report['transform'])
    vrt=ET.Element('VRTDataset',rasterXSize=str(report['width']),rasterYSize=str(report['height']))
    ET.SubElement(vrt,'SRS',dataAxisToSRSAxisMapping='2,1').text='EPSG:4326'
    ET.SubElement(vrt,'GeoTransform').text=','.join(map(str,transform.to_gdal()))
    band=ET.SubElement(vrt,'VRTRasterBand',dataType='Byte',band='1')
    ET.SubElement(band,'NoDataValue').text='255'
    for sheet in selected:
        if sheet['mask'] is None:
            continue
        source=root/'inputs/glo30-wbm'/sheet['mask']
        with rasterio.open(source) as src:
            expected=rasterio.Affine(*sheet['transform'])
            if (src.crs!=rasterio.crs.CRS.from_epsg(4326) or src.width!=sheet['width'] or src.height!=sheet['height']
                    or not np.allclose(src.transform,expected,atol=1e-10,rtol=0)):
                raise ValueError(f'River geometry requires a native-grid mask: {source.name}')
            x,y=(~transform)*(src.transform.c,src.transform.f)
            x,y=round(x),round(y)
            element=ET.SubElement(band,'ComplexSource')
            ET.SubElement(element,'SourceFilename',relativeToVRT='1').text=os.path.relpath(source,path.parent)
            ET.SubElement(element,'SourceBand').text='1'
            ET.SubElement(element,'NODATA').text='255'
            options=ET.SubElement(element,'OpenOptions')
            ET.SubElement(options,'OOI',key='OVERVIEW_LEVEL').text='NONE'
            ET.SubElement(element,'SrcRect',xOff='0',yOff='0',xSize=str(src.width),ySize=str(src.height))
            ET.SubElement(element,'DstRect',xOff=str(x),yOff=str(y),xSize=str(src.width),ySize=str(src.height))
    ET.indent(vrt)
    path.write_bytes(ET.tostring(vrt,encoding='utf-8',xml_declaration=True))


def geometry(root,label='final',names=None,halo=512):
    work=root/'rivers'/label
    (work/'geometry').mkdir(parents=True,exist_ok=True)
    inventory=read_json(root/'reports/inputs.json')
    selected=[s for s in inventory['sheets'] if names is None or s['name'] in names]
    if not selected or (names and len(selected)!=len(set(names))):
        raise ValueError('Missing or duplicate requested river input sheets')
    sources={}
    masks={m['name']:m['sha256'] for m in inventory['masks']}
    for sheet in selected:
        report=read_json(root/'reports/sheets'/(Path(sheet['name']).stem+'.json'))
        path=root/'sheets/coastal'/sheet['name']
        if sha256(path)!=report['output_sha256']:
            raise ValueError('Prepared elevation checksum differs')
        sources[sheet['name']]=report['output_sha256']
        if sheet['mask'] and sha256(root/'inputs/glo30-wbm'/sheet['mask'])!=masks[sheet['mask']]:
            raise ValueError('River classification input checksum differs')
    mask_recipe=dict(source='Copernicus GLO-30 WBM',river_class=3,ocean_class=1,
                     excluded_lake_class=2,
                     masks={s['mask']:masks[s['mask']] for s in selected if s['mask']})
    recipe=dict(sheets=sources,inputs_sha256=sha256(root/'reports/inputs.json'),halo=halo,
                river_source='glo30',river_mask=mask_recipe,
                allocation=allocation_recipe(root),
                geometry_sha256=sha256(Path(__file__).with_name('geometry.py')),
                staging_sha256=sha256(Path(__file__).with_name('staging.py')),
                driver_sha256=digest_json(inspect.getsource(geometry)+inspect.getsource(mask_mosaic)+inspect.getsource(atomic_npz)),
                numpy=np.__version__,gdal=rasterio.__gdal_version__,
                scipy=importlib.metadata.version('scipy'),scikit_image=importlib.metadata.version('scikit-image'))
    recipe_path=work/'geometry-recipe.json'
    if recipe_path.exists() and read_json(recipe_path)!=recipe:
        raise ValueError('Geometry inputs/code changed; use a different --label')
    atomic_json(recipe_path,recipe)
    height_path=work/'coastal.vrt';class_path=work/'classes.vrt'
    grid=create_mosaic([root/'sheets/coastal'/s['name'] for s in selected],height_path,work/'validity.vrt')
    mask_mosaic(root,selected,grid,class_path)
    atomic_json(work/'grid.json',grid)
    completed=[]
    with rasterio.Env(GDAL_CACHEMAX=64*1024**2,GDAL_NUM_THREADS='1',GDAL_MAX_DATASET_POOL_SIZE=24):
        with rasterio.open(height_path) as heights,rasterio.open(class_path) as classes:
            for i,sheet in enumerate(selected):
                started=time.monotonic()
                name=sheet['name']
                cache=work/'geometry'/(Path(name).stem+'.npz')
                record=cache.with_suffix('.json')
                if cache.exists():
                    report=read_json(record)
                    if report['sha256']!=sha256(cache) or report['recipe_sha256']!=digest_json(recipe):
                        raise ValueError('Geometry cache changed')
                elif len(selected)>10:
                    subprocess.run([sys.executable,'-m','river_flattening.staging',
                                    '--workdir',str(root),'--label',label,'--sheet',name],check=True)
                    report=read_json(record)
                    if report['sha256']!=sha256(cache) or report['recipe_sha256']!=digest_json(recipe):
                        raise ValueError('Geometry worker output differs')
                else:
                    transform=rasterio.Affine(*sheet['transform'])
                    x,y=(~heights.transform)*(transform.c,transform.f)
                    x,y=round(x),round(y)
                    x0=max(0,x-halo);y0=max(0,y-halo)
                    x1=min(heights.width,x+sheet['width']+halo);y1=min(heights.height,y+sheet['height']+halo)
                    window=Window(x0,y0,x1-x0,y1-y0)
                    h=heights.read(1,window=window)
                    mask=classes.read(1,window=window)
                    core=(y-y0,y-y0+sheet['height'],x-x0,x-x0+sheet['width'])
                    r0,r1,c0,c1=core
                    if np.any((h[r0:r1,c0:c1]!=-255)&(mask[r0:r1,c0:c1]>3)):
                        raise ValueError('Missing class coverage at a valid elevation')
                    lat=heights.transform.f+heights.transform.e*(y0+(y1-y0)/2)
                    dy=abs(heights.transform.e)*111195.0802
                    dx=heights.transform.a*111195.0802*math.cos(math.radians(lat))
                    data=extract_window(h,mask,dx=dx,dy=dy,grid_width=heights.width,
                                        row_offset=y0,col_offset=x0,core=core,
                                        allocator=make_allocator(root,work,name,recipe))
                    atomic_npz(cache,**data)
                    report=dict(sheet=name,sha256=sha256(cache),recipe_sha256=digest_json(recipe),
                                river_pixels=len(data['pixels']),skeleton_pixels=len(data['nodes']),
                                bank_samples=len(data['bank_heights']),elapsed_seconds=round(time.monotonic()-started,2))
                    atomic_json(record,report)
                    del h,mask,data
                completed.append(report)
                emit('river_geometry',complete=i+1,total=len(selected),**report)
                atomic_json(work/'geometry-status.json',dict(complete=i+1,total=len(selected),last_sheet=name))
    atomic_json(work/'geometry.json',dict(recipe=recipe,sheets=completed))
    return work


def main():
    import fcntl
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--label',default='final')
    parser.add_argument('--sheets',nargs='+')
    args=parser.parse_args()
    if Path(args.label).name!=args.label or args.label in ('.','..'):
        parser.error('Use a plain output label')
    root=setup(args.workdir)
    with (root/'logs/river-build.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        geometry(root,args.label,args.sheets)


if __name__=='__main__':main()
