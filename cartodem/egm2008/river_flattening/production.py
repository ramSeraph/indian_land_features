"""Assemble river geometry across sheets, fit once, then export per-sheet COGs.

Stages are restartable and keep input/output hashes. Parameters remain subject to
the method validation documented in this package; no external publication occurs.
"""
import argparse
import importlib.metadata
import inspect
import math
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

from build import atomic_json,create_mosaic,digest_json,emit,read_json,sha256,setup
from cog import export_cog
from .geometry import extract_window,compress_network
from .allocation import backend_recipe,make_allocator
from .study import fit_geometry
from .directions import fit_directed
from .joint import fit_joint
from .routing import route,sample_owners
from .boundaries import contact_files


def atomic_npz(path,**arrays):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.part')
    with temp.open('wb') as stream:
        np.savez_compressed(stream,**arrays)
    temp.replace(path)


def mask_mosaic(root,selected,report,path,mask_paths=None):
    transform=rasterio.Affine(*report['transform'])
    vrt=ET.Element('VRTDataset',rasterXSize=str(report['width']),rasterYSize=str(report['height']))
    ET.SubElement(vrt,'SRS',dataAxisToSRSAxisMapping='2,1').text='EPSG:4326'
    ET.SubElement(vrt,'GeoTransform').text=','.join(map(str,transform.to_gdal()))
    band=ET.SubElement(vrt,'VRTRasterBand',dataType='Byte',band='1')
    ET.SubElement(band,'NoDataValue').text='255'
    for sheet in selected:
        if sheet['mask'] is None:
            continue
        source=(mask_paths[sheet['name']] if mask_paths is not None else root/'inputs/glo30-wbm'/sheet['mask'])
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


def geometry(root,label='final',names=None,halo=512,river_source='glo30',allocation_backend='whitebox'):
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
    if river_source not in ('wris','glo30'):
        raise ValueError('Unknown river mask source')
    mask_paths=None
    mask_recipe=dict(source='Copernicus GLO-30 WBM',river_class=3,ocean_class=1,
                     excluded_lake_class=2,
                     masks={s['mask']:masks[s['mask']] for s in selected if s['mask']})
    if river_source=='wris':
        if len(selected)>10:
            # Reclaim GEOS/GDAL allocations between sheets. Reusing verified
            # masks below then assembles their ordinary shared manifest.
            from .staging import masks_in_processes
            masks_in_processes(root,label,selected)
        from .wris import prepare_masks
        mask_paths,mask_recipe=prepare_masks(root,work,selected)
    recipe=dict(sheets=sources,inputs_sha256=sha256(root/'reports/inputs.json'),halo=halo,
                river_source=river_source,river_mask=mask_recipe,
                allocation=backend_recipe(root,allocation_backend),
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
    mask_mosaic(root,selected,grid,class_path,mask_paths)
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
                    subprocess.run([sys.executable,'-m','river_flattening.staging','geometry',
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


def fit(root,label='final',spacing=180.,quantile=.2,stiffness=2.,max_support_gap=1200.,use_routing=True,join_method='joint'):
    work=root/'rivers'/label
    grid=read_json(work/'grid.json')
    if grid['width']*grid['height']>40_000_000:
        if join_method!='joint':raise ValueError('National fitting supports the agreed joint method only')
        from .national import fit as fit_national
        return fit_national(root,label,spacing,quantile,stiffness,max_support_gap,use_routing)
    records=read_json(work/'geometry.json')
    source_hashes={r['sheet']:r['sha256'] for r in records['sheets']}
    recipe=dict(geometry=source_hashes,spacing_m=spacing,quantile=quantile,stiffness=stiffness,
                join_method=join_method,
                maximum_support_gap_m=max_support_gap,
                direction_method='RivGraph guesses + robust CartoDEM bank trend + audited continuity',
                routing=read_json(work/'routing/report.json') if use_routing else None,
                code={p:sha256(Path(__file__).parent/p) for p in ('geometry.py','profiles.py','study.py','production.py','directions.py','joint.py','routing.py','boundaries.py','vendor/rivgraph.py','vendor/provenance.json')},
                versions={p:importlib.metadata.version(p) for p in ('numpy','scipy','clarabel','networkx')})
    target=work/'profiles.npz'
    if target.exists():
        old=read_json(work/'fit.json')
        if old['recipe']!=recipe or old['sha256']!=sha256(target):
            raise ValueError('Profile settings/code changed; use a different --label')
        return old
    collected={k:[] for k in ('nodes','bank_owners','bank_heights','ocean_owners')}
    for record in records['sheets']:
        path=work/'geometry'/(Path(record['sheet']).stem+'.npz')
        if sha256(path)!=record['sha256']:
            raise ValueError('Geometry cache checksum differs')
        with np.load(path) as data:
            for k in collected:
                collected[k].append(data[k])
    native={k:np.concatenate(parts) for k,parts in collected.items()}
    del collected
    native['nodes']=np.unique(native['nodes'])
    if not len(native['nodes']):
        raise ValueError('No river geometry in selected sheets')
    grid=read_json(work/'grid.json')
    transform=rasterio.Affine(*grid['transform'])
    # Local spherical ground distances vary by node latitude. The original
    # geographic raster grid is unchanged; these distances only set spacing and
    # regularization along the network.
    latitudes=transform.f+transform.e*(native['nodes']//grid['width']+.5)
    dy=abs(transform.e)*111195.0802
    dx=transform.a*111195.0802*np.cos(np.deg2rad(latitudes))
    network=compress_network(native['nodes'],grid['width'],dx,dy,spacing=spacing)
    emit('river_network',native_nodes=len(native['nodes']),stations=len(network['station_pixels']),
         reaches=len(network['reaches']))
    # Halo skeleton owners must be actual nodes of the assembled network.
    for key in ('bank_owners','ocean_owners'):
        j=np.searchsorted(network['ids'],native[key])
        okay=(j<len(network['ids']))&(network['ids'][np.minimum(j,len(network['ids'])-1)]==native[key])
        if not okay.all():
            raise ValueError(f'Halo geometry does not stitch: {int((~okay).sum())} missing {key}; increase halo and validate')
    routing=sample_owners(work,native,grid['width']) if use_routing else None
    if join_method not in ('joint','source-pins'):
        raise ValueError('Unknown river join method')
    fitter=fit_joint if join_method=='joint' else fit_directed
    network,components=fitter(native,network,spacing=spacing,quantile=quantile,stiffness=stiffness,
                                    max_support_gap=max_support_gap,routing=routing,width=grid['width'],
                                    contacts=contact_files(root,work))
    boundary=network['boundary_acceptance']
    atomic_json(work/'directions.json',dict(reaches=network['direction_diagnostics']))
    fields=('ids','left','right','fraction','station_pixels','fitted','bank_reference','bank_count','component','native_supported','supported_edges')
    atomic_npz(target,**{k:network[k] for k in fields})
    report=dict(recipe=recipe,sha256=sha256(target),components=components,
                native_nodes=len(network['ids']),stations=len(network['station_pixels']),
                fitted_stations=int(np.isfinite(network['fitted']).sum()),
                unsupported_skeleton_pixels=network['unsupported_skeleton_pixels'],
                direction_sha256=sha256(work/'directions.json'),boundary_acceptance=boundary,
                modeled_reaches=sum(d['profile_applied'] for d in network['direction_diagnostics']),
                joint_direction_hypotheses=sum(d['profile_applied'] and d.get('direction_source')=='joint profile hypothesis'
                                              for d in network['direction_diagnostics']),
                accepted_reaches=sum(d['status']=='accepted' for d in network['direction_diagnostics']),
                unresolved_reaches=sum(d['status']=='unresolved' for d in network['direction_diagnostics']))
    atomic_json(work/'fit.json',report)
    emit('river_fit',stations=report['stations'],fitted_stations=report['fitted_stations'],
         components=len(components),fitted_components=sum(r['status']=='fitted' for r in components))
    return report


def apply(root,label='final'):
    work=root/'rivers'/label
    geometry_report=read_json(work/'geometry.json')
    fitted=read_json(work/'fit.json')
    if fitted.get('partitioned'):
        from .national import apply as apply_national
        return apply_national(root,label)
    if fitted['sha256']!=sha256(work/'profiles.npz'):
        raise ValueError('Fitted river profiles changed')
    with np.load(work/'profiles.npz') as data:
        network={k:data[k] for k in data.files}
    grid=read_json(work/'grid.json')
    transform=rasterio.Affine(*grid['transform'])
    outdir=root/'output'/label/'cogs'
    outdir.mkdir(parents=True,exist_ok=True)
    reports=[]
    recipe=digest_json(dict(fit=fitted['sha256'],code=sha256(Path(__file__)),cog_code=sha256(Path(__file__).parents[1]/'cog.py')))
    for i,record in enumerate(geometry_report['sheets']):
        name=record['sheet']; output=outdir/name
        report_path=work/(Path(name).stem+'-cog.json')
        source=root/'sheets/coastal'/name
        if sha256(source)!=geometry_report['recipe']['sheets'][name]:
            raise ValueError('Coastal elevation input changed after river fitting')
        if output.exists():
            report=read_json(report_path)
            if report['recipe_sha256']!=recipe or report['sha256']!=sha256(output):
                raise ValueError('Existing final COG recipe or checksum differs')
            reports.append(report)
            continue
        with rasterio.open(source) as src:
            h=src.read(1)
            source_tags=src.tags(); profile=src.profile.copy()
            x,y=(~transform)*(src.transform.c,src.transform.f)
            x,y=round(x),round(y)
        path=work/'geometry'/(Path(name).stem+'.npz')
        if sha256(path)!=record['sha256']:
            raise ValueError('Geometry cache changed before application')
        with np.load(path) as data:
            pixels=data['pixels'];owners=data['owners']
        j=np.searchsorted(network['ids'],owners)
        if np.any(j>=len(network['ids'])) or np.any(network['ids'][j]!=owners):
            raise ValueError('A river pixel refers to a missing shared skeleton node')
        left,right=network['left'][j],network['right'][j]
        good=(left>=0)&(right>=0)&network['native_supported'][j]
        candidate=np.full(len(j),np.nan,dtype=np.float32)
        t=network['fraction'][j]
        candidate[good]=(1-t[good])*network['fitted'][left[good]]+t[good]*network['fitted'][right[good]]
        row=pixels//grid['width']-y;col=pixels%grid['width']-x
        if np.any((row<0)|(row>=h.shape[0])|(col<0)|(col>=h.shape[1])):
            raise ValueError('River geometry escaped its original sheet')
        original=h[row,col]
        valid=(original!=-255)&np.isfinite(original)
        treated=valid&np.isfinite(candidate)
        if np.any(candidate[treated]==-255):
            raise ValueError('River elevation collides with output nodata')
        delta=candidate[treated]-original[treated]
        nodata_before=int(np.sum(h==-255))
        h[row[treated],col[treated]]=candidate[treated]
        if int(np.sum(h==-255))!=nodata_before:
            raise ValueError('River processing changed original coverage')
        temporary=root/'tmp'/('river-'+name)
        try:
            with rasterio.open(temporary,'w',**profile) as dst:
                dst.write(h,1)
                dst.update_tags(**source_tags)
                dst.update_tags(RIVER_TREATMENT='Shared constrained CartoDEM bank profiles; unsupported components unchanged',
                                RIVER_TREATED_PIXELS=int(treated.sum()),RIVER_PROFILE_SHA256=fitted['sha256'],
                                RIVER_MASK=geometry_report['recipe']['river_source'],
                                NONRIVER_ELEVATIONS='Unchanged',VERTICAL_DATUM='EGM2008')
            report=export_cog(temporary,output)
        finally:
            temporary.unlink(missing_ok=True)
        report.update(recipe_sha256=recipe,source_sha256=geometry_report['recipe']['sheets'][name],
                      river_pixels=len(pixels),valid_river_pixels=int(valid.sum()),treated_river_pixels=int(treated.sum()),
                      nodata_pixels=nodata_before,
                      p95_absolute_river_change_m=float(np.percentile(np.abs(delta),95)) if len(delta) else None,
                      maximum_absolute_river_change_m=float(np.abs(delta).max()) if len(delta) else None)
        atomic_json(report_path,report)
        reports.append(report)
        emit('river_cog',complete=i+1,total=len(geometry_report['sheets']),**report)
    summary=dict(recipe_sha256=recipe,sheets=reports,nodata=-255,
                 river_source=geometry_report['recipe']['river_source'],
                 treated_river_pixels=sum(r['treated_river_pixels'] for r in reports),
                 valid_river_pixels=sum(r['valid_river_pixels'] for r in reports),
                 nodata_pixels=sum(r['nodata_pixels'] for r in reports))
    atomic_json(work/'cogs.json',summary)
    from .audit import audit
    audit(root,label)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['geometry','routing','fit','apply','run'])
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--label',default='final')
    parser.add_argument('--sheets',nargs='+')
    parser.add_argument('--halo',type=int,default=512)
    parser.add_argument('--allocation-backend',choices=['whitebox','scipy'],default='whitebox',
                        help='Raster cost allocation; SciPy is retained for small reference comparisons')
    parser.add_argument('--river-source',choices=['glo30','wris'],default='glo30',
                        help='GLO-30 class 3 covers the full DEM footprint; WRIS is retained for comparisons')
    parser.add_argument('--max-support-gap',type=float,default=1200.)
    parser.add_argument('--join-method',choices=['joint','source-pins'],default='joint',
                        help='Joint connected fit; source-pins retains the rejected per-pixel experiment')
    parser.add_argument('--bank-directions-only',action='store_true',help='Explicit ablation test without GRASS evidence')
    args=parser.parse_args()
    if Path(args.label).name!=args.label or args.halo<16:
        parser.error('Use a plain output label and a halo of at least 16 native pixels')
    root=setup(args.workdir)
    import fcntl
    with (root/'logs/river-build.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if args.stage in ('geometry','run'):
            geometry(root,args.label,args.sheets,args.halo,river_source=args.river_source,
                     allocation_backend=args.allocation_backend)
        if args.stage in ('routing','run') and not args.bank_directions_only:
            route(root,args.label)
        if args.stage in ('fit','run'):
            fit(root,args.label,max_support_gap=args.max_support_gap,use_routing=not args.bank_directions_only,
                join_method=args.join_method)
        if args.stage in ('apply','run'):
            apply(root,args.label)


if __name__=='__main__':
    main()
