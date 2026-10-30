"""Constrain joins to retained river heights and audit native-pixel contacts.

Checks the four cardinal pixel neighbours, including pairs across sheet borders.
Revert whole shared profile groups, never isolated pixels of a fitted cross-section.
The historical 10 m rejection gate remains for reproducing earlier experiments.
Current joins pin the correction to zero at the treated side of each interface.
This preserves the original local gradient; it does not repair source DEM steps.
"""
from pathlib import Path
import numpy as np
import rasterio
from build import read_json


def adjacent_contacts(pixels,owners,heights,width,sheets=None):
    order=np.argsort(pixels);pixels=pixels[order];owners=owners[order];heights=heights[order]
    sheets=sheets[order] if sheets is not None else None
    pairs=[];values=[]
    for offset in (1,width):
        j=np.searchsorted(pixels,pixels+offset);safe=np.minimum(j,max(0,len(j)-1))
        if not len(j): continue
        good=(j<len(j))&(pixels[safe]==pixels+offset)&(owners!=owners[safe])
        if offset==1: good&=pixels%width<width-1
        good&=np.isfinite(heights)&np.isfinite(heights[safe])&(heights!=-255)&(heights[safe]!=-255)
        if sheets is not None: good&=sheets!=sheets[safe]
        i=np.flatnonzero(good);j=j[good]
        pairs.append(np.column_stack((owners[i],owners[j])))
        values.append(np.column_stack((heights[i],heights[j])))
    if not pairs: return dict(owners=np.empty((0,2),int),heights=np.empty((0,2),'float32'))
    return dict(owners=np.concatenate(pairs),heights=np.concatenate(values))


def contact_files(root,work):
    from .production import atomic_npz
    records=read_json(work/'geometry.json')['sheets'];grid=read_json(work/'grid.json')
    t=rasterio.Affine(*grid['transform']);width=grid['width'];paths=[];border=[]
    folder=work/'contacts';folder.mkdir(exist_ok=True)
    for k,record in enumerate(records):
        name=record['sheet']
        with np.load(work/'geometry'/(Path(name).stem+'.npz')) as data:
            pixels,owners=data['pixels'],data['owners']
        with rasterio.open(root/'sheets/coastal'/name) as src:
            x,y=(~t)*(src.transform.c,src.transform.f);x,y=round(x),round(y)
            row=pixels//width-y;col=pixels%width-x
            heights=src.read(1)[row,col]
            edge=(row==0)|(row==src.height-1)|(col==0)|(col==src.width-1)
        path=folder/(Path(name).stem+'.npz')
        atomic_npz(path,**adjacent_contacts(pixels,owners,heights,width));paths.append(path)
        border.append((pixels[edge],owners[edge],heights[edge],np.full(edge.sum(),k,int)))
    if border:
        pixels,owners,heights,sheets=(np.concatenate([b[i] for b in border]) for i in range(4))
        path=folder/'sheet-borders.npz'
        atomic_npz(path,**adjacent_contacts(pixels,owners,heights,width,sheets));paths.append(path)
    return paths


def native_values(network):
    left,right=network['left'],network['right'];t=network['fraction']
    good=(left>=0)&(right>=0)&network['native_supported']
    values=np.full(len(left),np.nan,dtype=np.float32)
    values[good]=(1-t[good])*network['fitted'][left[good]]+t[good]*network['fitted'][right[good]]
    groups=np.full(len(left),-1,int)
    groups[good]=network['component'][left[good]]
    return values,groups


def read_contacts(table):
    if isinstance(table,(str,Path)):
        with np.load(table) as data:
            return data['owners'],data['heights']
    return table['owners'],table['heights']


def join_samples(network,tables):
    """Anchor active/retained and independently fitted group interfaces.

    Samples refer to the treated pixel's own pre-fit height. Keeping its correction
    zero at the interface preserves the source's two-centre slope, instead of
    imposing a flat step or copying the neighbouring pixel's elevation.
    """
    values,groups=native_values(network);ids=network['ids'];samples={};pairs=0
    for table in tables:
        owners,original=read_contacts(table)
        if not len(owners):continue
        index=np.searchsorted(ids,owners)
        if np.any(index>=len(ids)) or np.any(ids[index]!=owners):
            raise ValueError('Boundary contact has a missing owner')
        active=np.isfinite(values[index]);group=groups[index]
        boundary=(active[:,0]!=active[:,1])|(active.all(axis=1)&(group[:,0]!=group[:,1]))
        pairs+=int(boundary.sum())
        for side in (0,1):
            take=boundary&active[:,side]
            local=index[take,side]
            rows=np.column_stack((network['left'][local],network['right'][local],
                                  network['fraction'][local],original[take,side]))
            for g in np.unique(group[take,side]):
                samples.setdefault(int(g),[]).append(rows[group[take,side]==g])
    # Exact duplicate pixel constraints add no information and can be numerous
    # along a wide channel. Distinct source heights remain distinct constraints.
    return {g:np.unique(np.concatenate(rows),axis=0) for g,rows in samples.items()},pairs


def constrain_joins(network,tables,refit):
    """Refit with boundary constraints; retain whole incompatible groups.

    refit(group, samples) updates that group's station values only after its
    constrained solver and support checks succeed. Dropping a group exposes new
    interfaces, so rescan and refit affected neighbours until the set is stable.
    """
    tables=list(tables);rejected={};iterations=[]
    for _ in range(len(network['fitted'])+1):
        samples,pairs=join_samples(network,tables)
        failures={}
        for group,rows in samples.items():
            try:
                refit(group,rows)
            except ValueError as exc:
                failures[group]=str(exc)
        iterations.append(dict(boundary_pairs=pairs,constrained_groups=len(samples),
                               constraints=sum(len(r) for r in samples.values()),
                               rejected_groups={str(k):v for k,v in failures.items()}))
        if not failures:
            worst=0.
            for rows in samples.values():
                a,b=rows[:,:2].astype(int).T;t=rows[:,2]
                predicted=((1-t)*network['fitted'][a]+t*network['fitted'][b]).astype(np.float32)
                worst=max(worst,float(np.abs(predicted.astype(float)-rows[:,3]).max(initial=0)))
            if worst>.001:
                raise ValueError('Stored Float32 profiles do not preserve treatment-boundary levels')
            return dict(method='Boundary-constrained refit; zero correction at retained interfaces',
                        numerical_tolerance_m=.001,maximum_boundary_mismatch_m=worst,
                        boundary_pairs=pairs,constrained_groups=len(samples),
                        rejected_profile_groups=sorted(rejected),
                        rejection_reasons={str(k):v for k,v in rejected.items()},iterations=iterations)
        rejected.update(failures)
        drop=np.isin(network['component'],list(failures))
        network['fitted'][drop]=np.nan
        left,right=network['left'],network['right'];good=(left>=0)&(right>=0)
        network['native_supported'][good]&=~(drop[left[good]]|drop[right[good]])
        edges=network['supported_edges']
        network['supported_edges']=edges[~(drop[edges[:,0]]|drop[edges[:,1]])]
    raise ValueError('River boundary constraints did not converge')


def guard(network,tables,maximum_added_step=10.):
    """Tables may be dictionaries (tests) or NPZ paths (bounded-memory build)."""
    ids=network['ids'];left=network['left'];right=network['right'];t=network['fraction']
    iterations=[];rejected=set()
    for _ in range(len(network['fitted'])+1):
        good=(left>=0)&(right>=0)&network['native_supported']
        candidate=np.full(len(ids),np.nan,dtype=np.float32)
        candidate[good]=(1-t[good])*network['fitted'][left[good]]+t[good]*network['fitted'][right[good]]
        groups=np.full(len(ids),-1,int);groups[good]=network['component'][left[good]]
        bad_groups=set();bad_pairs=0;worst=0.
        for table in tables:
            if isinstance(table,(str,Path)):
                with np.load(table) as data: owners,original=data['owners'],data['heights']
            else: owners,original=table['owners'],table['heights']
            if not len(owners): continue
            index=np.searchsorted(ids,owners)
            if np.any(index>=len(ids)) or np.any(ids[index]!=owners): raise ValueError('Boundary contact has a missing owner')
            fitted=candidate[index];applied=np.isfinite(fitted)
            after=np.where(applied,fitted,original)
            before_step=np.abs(original[:,1]-original[:,0]);after_step=np.abs(after[:,1]-after[:,0])
            increase=after_step-before_step
            bad=(after_step>maximum_added_step)&(increase>maximum_added_step)
            bad_pairs+=int(bad.sum());worst=max(worst,float(np.max(increase,initial=0)))
            changed=applied&(after!=original)&bad[:,None]
            bad_groups.update(int(g) for g in groups[index][changed] if g>=0)
        if not bad_groups:
            return dict(maximum_added_step_m=maximum_added_step,rejected_profile_groups=sorted(rejected),
                        iterations=iterations,remaining_bad_pairs=bad_pairs,maximum_remaining_added_step_m=worst)
        iterations.append(dict(rejected_groups=sorted(bad_groups),bad_pairs=bad_pairs,maximum_added_step_m=worst))
        rejected.update(bad_groups)
        drop=np.isin(network['component'],list(bad_groups))
        network['fitted'][drop]=np.nan
        network['native_supported'][good]&=~(drop[left[good]]|drop[right[good]])
        edges=network['supported_edges'];network['supported_edges']=edges[~(drop[edges[:,0]]|drop[edges[:,1]])]
    raise ValueError('River boundary acceptance gate did not converge')
