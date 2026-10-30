"""Verify a river pilot independently against its coastal baseline and masks."""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window,from_bounds

from build import atomic_json,read_json,sha256


def audit_native_joins(pixels,before,after,applied,groups,width,sheets=None):
    """Inspect written raster values at treatment interfaces, including seams."""
    order=np.argsort(pixels)
    pixels,before,after,applied,groups=(v[order] for v in (pixels,before,after,applied,groups))
    if sheets is not None:sheets=sheets[order]
    pairs=0;worst=0.;maximum_signed_step_change=0.
    for offset in (1,width):
        if not len(pixels):continue
        j=np.searchsorted(pixels,pixels+offset);safe=np.minimum(j,len(j)-1)
        good=(j<len(j))&(pixels[safe]==pixels+offset)
        if offset==1:good&=pixels%width<width-1
        if sheets is not None:good&=sheets!=sheets[safe]
        good&=(before!=-255)&(before[safe]!=-255)&np.isfinite(before)&np.isfinite(before[safe])
        boundary=(applied!=applied[safe])|(applied&applied[safe]&(groups!=groups[safe]))
        i=np.flatnonzero(good&boundary);jj=j[i]
        pairs+=len(i)
        for side in (i,jj):
            changed=side[applied[side]]
            worst=max(worst,float(np.abs(after[changed].astype(float)-before[changed]).max(initial=0)))
        step_change=(after[jj].astype(float)-after[i])-(before[jj].astype(float)-before[i])
        maximum_signed_step_change=max(maximum_signed_step_change,float(np.abs(step_change).max(initial=0)))
    return dict(boundary_pairs=pairs,maximum_treated_boundary_correction_m=worst,
                maximum_signed_step_change_m=maximum_signed_step_change)


def audit(root,label,previous=None):
    root=Path(root);work=root/'rivers'/label
    geometry=read_json(work/'geometry.json');fit=read_json(work/'fit.json')
    grid=read_json(work/'grid.json');transform=rasterio.Affine(*grid['transform'])
    with np.load(work/'profiles.npz') as data:
        values=data['fitted'];edges=data['supported_edges'];ids=data['ids']
        left,right,supported=data['left'],data['right'],data['native_supported']
        labels=data['component']
    active=(left>=0)&(right>=0)&supported
    active[active]&=np.isfinite(values[left[active]])&np.isfinite(values[right[active]])
    groups=np.full(len(ids),-1,int);groups[active]=labels[left[active]]
    if len(edges) and np.any(values[edges[:,0]]<values[edges[:,1]]):
        raise ValueError('A stored Float32 profile rises downstream')
    old=None
    if previous:
        with np.load(root/'rivers'/previous/'profiles.npz') as data:
            if not np.array_equal(ids,data['ids']): raise ValueError('Comparison skeleton differs')
            station=np.where(data['fraction']<=.5,data['left'],data['right'])
            old=np.where(station>=0,data['component'][np.maximum(station,0)],-1)
    component_counts={};sheets=[];borders=[];joins=[]
    with rasterio.open(work/'classes.vrt') as classes:
        for item in geometry['sheets']:
            name=item['sheet'];source=root/'sheets/coastal'/name;output=root/'output'/label/'cogs'/name
            with rasterio.open(source) as src,rasterio.open(output) as dst:
                before=src.read(1);after=dst.read(1)
                assert src.transform==dst.transform and src.shape==dst.shape and src.crs==dst.crs
                x,y=(~transform)*(src.transform.c,src.transform.f);x,y=round(x),round(y)
                mask=classes.read(1,window=Window(x,y,src.width,src.height))
                valid=(before!=-255)&np.isfinite(before);river=valid&(mask==3)
                preserve=~river
                if not np.array_equal(before[preserve],after[preserve],equal_nan=True):
                    raise ValueError('A non-river, ocean, lake or nodata pixel changed')
                if not np.array_equal(src.read_masks(1),dst.read_masks(1)):
                    raise ValueError('Native validity changed')
                delta=after-before;changed=river&(delta!=0)
                ii=np.flatnonzero(changed);take=ii[np.argsort(np.abs(delta.ravel()[ii]))[-20:][::-1]]
                largest=[]
                for p in take:
                    lon,lat=src.transform*(int(p%src.width)+.5,int(p//src.width)+.5)
                    largest.append(dict(latitude=lat,longitude=lon,before_m=float(before.ravel()[p]),after_m=float(after.ravel()[p])))
                sharp=[]
                for dr,dc in [(0,1),(1,0)]:
                    a=np.s_[:src.height-dr,:src.width-dc];b=np.s_[dr:,dc:]
                    water=river[a]&river[b]
                    jump0=np.abs(before[b]-before[a]);jump1=np.abs(after[b]-after[a])
                    bad=water&(jump1>10)&(jump1-jump0>10)
                    locations=[]
                    for row,col in np.argwhere(bad):
                        lon,lat=src.transform*(int(col)+.5,int(row)+.5)
                        locations.append(dict(latitude=lat,longitude=lon,row=int(row),col=int(col),
                                              before=[float(before[row,col]),float(before[row+dr,col+dc])],
                                              after=[float(after[row,col]),float(after[row+dr,col+dc])]))
                    sharp.append(dict(offset=[dr,dc],new_jump_over_10m=int(np.sum(water&(jump1>10)&(jump1-jump0>10))),
                                      maximum_added_jump_m=float(np.max(np.where(water,jump1-jump0,0))),locations=locations))
                sheets.append(dict(sheet=name,output_sha256=sha256(output),valid_river_pixels=int(river.sum()),
                                   changed_river_pixels=int(changed.sum()),native_validity_unchanged=True,
                                   nonriver_and_nodata_unchanged=True,largest_changes=largest,adjacent_water_jumps=sharp))
                with np.load(work/'geometry'/(Path(name).stem+'.npz')) as data:
                    pixels,owners=data['pixels'],data['owners']
                row=pixels//grid['width']-y;col=pixels%grid['width']-x
                node=np.searchsorted(ids,owners)
                treated=active[node]&valid[row,col]
                original,actual=before[row,col],after[row,col]
                if not np.array_equal(original[~treated],actual[~treated],equal_nan=True):
                    raise ValueError('An unresolved river pixel changed')
                local_groups=groups[node]
                check=audit_native_joins(pixels,original,actual,treated,local_groups,grid['width'])
                joins.append(dict(sheet=name,**check))
                border=(row==0)|(row==src.height-1)|(col==0)|(col==src.width-1)
                borders.append(tuple(v[border] for v in (pixels,original,actual,treated,local_groups))+
                               (np.full(int(border.sum()),len(sheets)-1,int),))
                if old is not None:
                    group=old[np.searchsorted(ids,owners)]
                    for k in np.unique(group):
                        c=component_counts.setdefault(str(k),dict(valid=0,changed=0))
                        use=group==k;c['valid']+=int(valid[row[use],col[use]].sum());c['changed']+=int(changed[row[use],col[use]].sum())
    if borders:
        p,b,a,t,g,s=(np.concatenate([row[i] for row in borders]) for i in range(6))
        joins.append(dict(sheet='Across sheet borders',**audit_native_joins(p,b,a,t,g,grid['width'],s)))
    if fit['boundary_acceptance'].get('method','').startswith('Boundary-constrained'):
        if any(r['maximum_treated_boundary_correction_m']>.001 for r in joins):
            raise ValueError('Written COG has an unconstrained treatment boundary')
    if fit['boundary_acceptance'].get('method','').startswith('Shared junction'):
        if any(r['boundary_pairs'] for r in joins):
            raise ValueError('Written COG contains a partial-treatment interface inside a connected river')
    directions=read_json(work/'directions.json')['reaches']
    report=dict(label=label,auditor_sha256=sha256(Path(__file__)),sheets=sheets,downstream_profile_violations=0,
                accepted_reaches=fit['accepted_reaches'],unresolved_reaches=fit['unresolved_reaches'],
                evidence=dict(Counter(e for r in directions if r['status']=='accepted' for e in r['evidence'])),
                direction_rejections=dict(Counter(r.get('reason','Insufficient direction evidence') for r in directions if r['status']=='unresolved')),
                fit_rejections=dict(Counter(r['reason'] for r in fit['components'] if r['status']!='fitted')),
                native_treatment_joins=joins,
                previous_component_counts=component_counts)
    atomic_json(work/'audit.json',report)
    return report


def crop(root,label,bounds,path):
    work=Path(root)/'rivers'/label
    with rasterio.open(work/'coastal.vrt') as src,rasterio.open(work/'classes.vrt') as mask:
        w=from_bounds(*bounds,transform=src.transform).round_offsets().round_lengths()
        baseline=src.read(1,window=w);classes=mask.read(1,window=w);t=src.window_transform(w)
    # The diagnostic crop must fit one output sheet, checked explicitly.
    for p in (Path(root)/'output'/label/'cogs').glob('*.tif'):
        with rasterio.open(p) as src:
            if src.bounds.left<=bounds[0] and src.bounds.right>=bounds[2] and src.bounds.bottom<=bounds[1] and src.bounds.top>=bounds[3]:
                window=from_bounds(*bounds,transform=src.transform).round_offsets().round_lengths()
                after=src.read(1,window=window)
                np.savez_compressed(path,before=baseline,after=after,classes=classes,transform=np.array(tuple(t)[:6]))
                return
    raise ValueError('Diagnostic crop does not fit one final sheet')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,required=True);parser.add_argument('--label',required=True)
    parser.add_argument('--previous');parser.add_argument('--crop',type=float,nargs=4);parser.add_argument('--crop-output',type=Path)
    args=parser.parse_args()
    print(audit(args.workdir,args.label,args.previous))
    if args.crop: crop(args.workdir,args.label,args.crop,args.crop_output)
