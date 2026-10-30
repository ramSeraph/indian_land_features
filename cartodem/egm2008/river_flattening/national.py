"""Fit complete connected river components with disk-backed staging.

Processing boundaries are not river boundaries. A sparse national skeleton joins
all sheets, including valid native water contacts, before components are staged
in SQLite and fitted separately. Each component uses the existing joint solver.
"""
from collections import Counter
import io
import json
import math
from pathlib import Path
import sqlite3
import time

import numpy as np
import rasterio
from rasterio.windows import Window
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from build import atomic_json,digest_json,emit,read_json,sha256
from cog import export_cog
from .geometry import compress_network
from .joint import fit_joint
from .boundaries import contact_files,native_values
from validation.river_audit import audit_native_joins


def pack(data):
    stream=io.BytesIO();np.savez_compressed(stream,**data)
    return stream.getvalue()


def unpack(blob):
    with np.load(io.BytesIO(blob)) as data:return {k:data[k] for k in data.files}


def lookup(ids,query):
    index=np.searchsorted(ids,query)
    if np.any(index>=len(ids)) or np.any(ids[index]!=query):
        raise ValueError('Halo geometry does not stitch: an owner is absent from the national skeleton')
    return index


def skeleton_components(ids,width):
    """Boolean CSR connectivity avoids retaining metric graphs for all India."""
    n=len(ids)
    if not n:return 0,np.empty(0,np.int32)
    rows=[];cols=[]
    for dr,dc in ((0,1),(1,-1),(1,0),(1,1)):
        target=ids+dr*width+dc;index=np.searchsorted(ids,target)
        good=(index<n)&(ids[np.minimum(index,n-1)]==target)
        good&=(ids%width+dc>=0)&(ids%width+dc<width)
        a=np.flatnonzero(good).astype('int32');b=index[good].astype('int32')
        rows.extend((a,b));cols.extend((b,a))
    r=np.concatenate(rows);c=np.concatenate(cols)
    graph=sparse.csr_matrix((np.ones(len(r),np.uint8),(r,c)),shape=(n,n))
    return connected_components(graph,directed=False)


def root_of(parent,a):
    while parent[a]!=a:
        parent[a]=parent[parent[a]];a=int(parent[a])
    return a


def pieces(groups,arrays,result):
    """Group one sheet's bounded arrays; equal source samples occur only once."""
    if not len(groups):return
    order=np.argsort(groups,kind='stable');values=groups[order]
    starts=np.r_[0,np.flatnonzero(np.diff(values))+1,len(values)]
    for a,b in zip(starts[:-1],starts[1:]):
        item=result.setdefault(int(values[a]),{})
        take=order[a:b]
        for name,data in arrays.items():item[name]=data[take]


def partition(root,label,use_routing=True):
    root=Path(root);work=root/'rivers'/label;folder=work/'components';folder.mkdir(exist_ok=True)
    records=read_json(work/'geometry.json');grid=read_json(work/'grid.json')
    routing=read_json(work/'routing/report.json') if use_routing else None
    if routing is not None and routing.get('mode')!='regional':
        raise ValueError('Partitioned fitting requires regional direction evidence')
    recipe=dict(geometry_sha256=sha256(work/'geometry.json'),routing=routing,
                code_sha256=sha256(Path(__file__)),native_valid_water_contacts=True)
    done=folder/'partition.json';database=folder/'chunks.sqlite'
    if done.exists():
        report=read_json(done)
        if report['recipe']!=recipe:raise ValueError('Partition recipe changed; use a new label')
        for name,digest in report['files'].items():
            if sha256(folder/name)!=digest:raise ValueError('Partition data changed')
        return report
    parts=[]
    for record in records['sheets']:
        path=work/'geometry'/(Path(record['sheet']).stem+'.npz')
        if sha256(path)!=record['sha256']:raise ValueError('Geometry cache changed')
        with np.load(path) as data:parts.append(data['nodes'])
    ids=np.unique(np.concatenate(parts));del parts
    count,labels=skeleton_components(ids,grid['width'])
    parent=np.arange(count,dtype=np.int32)
    tables=contact_files(root,work)
    for path in tables:
        with np.load(path) as data:owners=data['owners']
        if not len(owners):continue
        pairs=np.unique(np.sort(labels[lookup(ids,owners)],axis=1),axis=0)
        for a,b in pairs:
            a=root_of(parent,int(a));b=root_of(parent,int(b))
            if a!=b:parent[max(a,b)]=min(a,b)
    for i in range(count):parent[i]=root_of(parent,i)
    labels=parent[labels]
    np.save(folder/'nodes.npy',ids);np.save(folder/'labels.npy',labels)
    pending=database.with_name(database.name+'.part');pending.unlink(missing_ok=True)
    regions={r['sheet']:r for r in routing['regions']} if routing else {}
    with sqlite3.connect(pending) as db:
        db.execute('PRAGMA journal_mode=DELETE');db.execute('PRAGMA cache_size=-32768')
        db.execute('CREATE TABLE chunks (component INTEGER, source TEXT, data BLOB, PRIMARY KEY(component,source))')
        db.execute('CREATE TABLE jobs (component INTEGER PRIMARY KEY, native_nodes INTEGER NOT NULL)')
        keys,counts=np.unique(labels,return_counts=True)
        db.executemany('INSERT INTO jobs VALUES (?,?)',((int(k),int(n)) for k,n in zip(keys,counts)))
        for number,record in enumerate(records['sheets']):
            name=record['sheet'];chunks={}
            with np.load(work/'geometry'/(Path(name).stem+'.npz')) as data:
                nodes=data['nodes'];bank=data['bank_owners'];ocean=data['ocean_owners']
                pieces(labels[lookup(ids,nodes)],dict(nodes=nodes),chunks)
                pieces(labels[lookup(ids,bank)],dict(bank_owners=bank,bank_heights=data['bank_heights']),chunks)
                pieces(labels[lookup(ids,ocean)],dict(ocean_owners=ocean),chunks)
            if routing:
                r=regions[name];path=work/'routing'/r['evidence']
                if sha256(path)!=r['sha256']:raise ValueError('Regional owner evidence changed')
                with np.load(path) as data:
                    nodes=data['nodes']
                    pieces(labels[lookup(ids,nodes)],{('route_'+k):data[k] for k in data.files},chunks)
            # Contacts are reduced to unique owner pairs. The joint fit uses
            # their topology; it never pins the raw pixel heights.
            with np.load(tables[number]) as data:pairs=np.unique(data['owners'],axis=0)
            if len(pairs):pieces(labels[lookup(ids,pairs[:,0])],dict(contacts=pairs),chunks)
            db.executemany('INSERT INTO chunks VALUES (?,?,?)',
                           ((component,name,pack(data)) for component,data in chunks.items()))
            db.commit()
            emit('river_partition',complete=number+1,total=len(records['sheets']),sheet=name,components=len(chunks))
        if len(tables)>len(records['sheets']):
            with np.load(tables[-1]) as data:pairs=np.unique(data['owners'],axis=0)
            chunks={}
            if len(pairs):pieces(labels[lookup(ids,pairs[:,0])],dict(contacts=pairs),chunks)
            db.executemany('INSERT INTO chunks VALUES (?,?,?)',((k,'sheet-borders',pack(v)) for k,v in chunks.items()))
            db.commit()
    pending.replace(database)
    report=dict(recipe=recipe,native_nodes=len(ids),components=len(keys),
                largest_component_nodes=int(counts.max(initial=0)),
                files={p:sha256(folder/p) for p in ('nodes.npy','labels.npy','chunks.sqlite')})
    atomic_json(done,report)
    return report


def fit(root,label,spacing=180.,quantile=.2,stiffness=2.,max_support_gap=1200.,use_routing=True):
    root=Path(root);work=root/'rivers'/label;folder=work/'components'
    staged=partition(root,label,use_routing=use_routing)
    recipe=dict(partition_sha256=sha256(folder/'partition.json'),spacing_m=spacing,quantile=quantile,
                stiffness=stiffness,maximum_support_gap_m=max_support_gap,
                method='existing joint profiles; complete global components',
                code={p:sha256(Path(__file__).with_name(p)) for p in
                      ('national.py','joint.py','profiles.py','geometry.py','directions.py','boundaries.py')})
    recipe_id=digest_json(recipe);target=folder/'models.sqlite';done=work/'fit.json'
    if done.exists():
        report=read_json(done)
        if report['recipe']!=recipe or report['sha256']!=sha256(target):raise ValueError('National fit recipe or models changed')
        return report
    grid=read_json(work/'grid.json');t=rasterio.Affine(*grid['transform'])
    with sqlite3.connect(f'file:{folder/"chunks.sqlite"}?mode=ro',uri=True) as source,sqlite3.connect(target) as db:
        db.execute('PRAGMA cache_size=-32768')
        db.execute('CREATE TABLE IF NOT EXISTS metadata (recipe TEXT NOT NULL)')
        old=db.execute('SELECT recipe FROM metadata').fetchone()
        if old and old[0]!=recipe_id:raise ValueError('Checkpointed national fit settings changed')
        if old is None:db.execute('INSERT INTO metadata VALUES (?)',(recipe_id,))
        db.execute('CREATE TABLE IF NOT EXISTS models (component INTEGER PRIMARY KEY, report TEXT NOT NULL, data BLOB)')
        db.commit()
        jobs=source.execute('SELECT component,native_nodes FROM jobs ORDER BY component').fetchall()
        complete={row[0] for row in db.execute('SELECT component FROM models')}
        for component,size in jobs:
            if component in complete:continue
            started=time.monotonic();collected={}
            for (blob,) in source.execute('SELECT data FROM chunks WHERE component=? ORDER BY source',(component,)):
                for key,array in unpack(blob).items():collected.setdefault(key,[]).append(array)
            data={k:np.concatenate(v) for k,v in collected.items()};del collected
            native={k:data.get(k,np.empty(0,np.float32 if k=='bank_heights' else np.int64)) for k in
                    ('nodes','bank_owners','bank_heights','ocean_owners')}
            native['nodes']=np.unique(native['nodes']);nodes=native['nodes']
            if len(nodes)!=size:raise ValueError('Staged component node count differs')
            lat=t.f+t.e*(nodes//grid['width']+.5)
            dx=t.a*111195.0802*np.cos(np.deg2rad(lat));dy=abs(t.e)*111195.0802
            network=compress_network(nodes,grid['width'],dx,dy,spacing=spacing)
            routing=None
            if use_routing:
                routing=dict(vector=np.zeros((len(nodes),2)),weight=np.zeros(len(nodes)),count=np.zeros(len(nodes),np.int64))
                if len(data.get('route_nodes',[])):
                    index=lookup(nodes,data['route_nodes'])
                    for axis in range(2):np.add.at(routing['vector'][:,axis],index,data['route_vector'][:,axis])
                    for key in ('weight','count'):np.add.at(routing[key],index,data['route_'+key])
            pairs=data.get('contacts',np.empty((0,2),np.int64))
            contacts=[dict(owners=pairs,heights=np.zeros(pairs.shape,np.float32))]
            # An isolated skeleton pixel has no reach, but remains an explicit
            # unsupported component and must not break the national run.
            if not len(network['reaches']):
                report=dict(components=[dict(status='unchanged',reason='No modelled river reach')],
                            native_nodes=len(nodes),stations=len(network['station_pixels']),
                            fitted_stations=0,accepted_reaches=0,unresolved_reaches=0,
                            modeled_reaches=0,joint_direction_hypotheses=0)
                model=None
            else:
                network,components=fit_joint(native,network,routing=routing,width=grid['width'],
                                             spacing=spacing,quantile=quantile,stiffness=stiffness,
                                             max_support_gap=max_support_gap,contacts=contacts)
                heights,groups=native_values(network)
                groups=np.where(groups>=0,np.int64(component)*2**32+groups,-1)
                directions=network['direction_diagnostics']
                report=dict(components=components,native_nodes=len(nodes),stations=len(network['station_pixels']),
                            fitted_stations=int(np.isfinite(network['fitted']).sum()),
                            accepted_reaches=sum(d['status']=='accepted' for d in directions),
                            unresolved_reaches=sum(d['status']=='unresolved' for d in directions),
                            modeled_reaches=sum(d['profile_applied'] for d in directions),
                            joint_direction_hypotheses=sum(d['profile_applied'] and d.get('direction_source')=='joint profile hypothesis' for d in directions),
                            directions=directions,boundary_acceptance=network['boundary_acceptance'])
                # Store station arrays for independent downstream validation,
                # and native projections for bounded per-sheet application.
                model=pack(dict(ids=nodes,heights=heights,groups=groups,fitted=network['fitted'],
                                supported_edges=network['supported_edges'],station_pixels=network['station_pixels'],
                                left=network['left'],right=network['right'],fraction=network['fraction']))
            report.update(component=component,elapsed_seconds=round(time.monotonic()-started,2))
            db.execute('INSERT INTO models VALUES (?,?,?)',(component,json.dumps(report),model));db.commit()
            complete.add(component)
            if len(complete)%100==0 or len(complete)==len(jobs) or report['elapsed_seconds']>10:
                emit('river_component_fit',complete=len(complete),total=len(jobs),component=component,
                     nodes=size,fitted_stations=report['fitted_stations'])
                atomic_json(folder/'fit-status.json',dict(complete=len(complete),total=len(jobs),component=component))
        totals=Counter();reasons=Counter()
        for (text,) in db.execute('SELECT report FROM models'):
            row=json.loads(text)
            for key in ('native_nodes','stations','fitted_stations','accepted_reaches','unresolved_reaches','modeled_reaches','joint_direction_hypotheses'):
                totals[key]+=row[key]
            for c in row['components']:
                if c['status']!='fitted':reasons[c['reason']]+=1
    report=dict(partitioned=True,recipe=recipe,sha256=sha256(target),components=len(jobs),
                fit_rejections=dict(reasons),**totals,
                boundary_acceptance=dict(method='Shared junction variables; complete connected treatment',raw_height_pins=False))
    atomic_json(done,report)
    return report


def project_owners(owners,ids,labels,db):
    groups=labels[lookup(ids,owners)]
    values=np.full(len(owners),np.nan,np.float32);profile_groups=np.full(len(owners),-1,np.int64)
    for group in np.unique(groups):
        row=db.execute('SELECT data FROM models WHERE component=?',(int(group),)).fetchone()
        if row is None:raise ValueError('A river component has no completed fit record')
        if row[0] is None:continue
        model=unpack(row[0]);edges=model['supported_edges'];fitted=model['fitted']
        if len(edges) and np.any(fitted[edges[:,0]]<fitted[edges[:,1]]):
            raise ValueError('Stored component profile rises downstream')
        selected=np.flatnonzero(groups==group);index=lookup(model['ids'],owners[selected])
        values[selected]=model['heights'][index];profile_groups[selected]=model['groups'][index]
    return values,profile_groups


def gradient_report(before,after,river,transform):
    reports=[]
    for dr,dc in ((0,1),(1,0)):
        a=np.s_[:before.shape[0]-dr,:before.shape[1]-dc];b=np.s_[dr:,dc:]
        water=river[a]&river[b]
        original=np.abs(before[b]-before[a]);actual=np.abs(after[b]-after[a])
        increase=np.where(water,actual-original,0)
        flat=increase.ravel();n=min(10,len(flat))
        selected=np.argpartition(flat,-n)[-n:] if n else []
        selected=sorted((int(i) for i in selected if flat[i]>0),key=lambda i:float(flat[i]),reverse=True)
        examples=[]
        for i in selected:
            row,col=divmod(i,increase.shape[1]);lon,lat=transform*(col+.5,row+.5)
            examples.append(dict(latitude=lat,longitude=lon,added_jump_m=float(flat[i]),
                                 before=[float(before[row,col]),float(before[row+dr,col+dc])],
                                 after=[float(after[row,col]),float(after[row+dr,col+dc])]))
        reports.append(dict(offset=[dr,dc],maximum_added_jump_m=float(increase.max(initial=0)),
                            maximum_water_jump_m=float(np.where(water,actual,0).max(initial=0)),
                            new_jump_over_10m=int(np.sum(water&(actual>10)&(increase>10))),
                            largest_added_jumps=examples))
    return reports


def apply(root,label):
    root=Path(root);work=root/'rivers'/label;folder=work/'components'
    geometry=read_json(work/'geometry.json');fitted=read_json(work/'fit.json')
    partitioned=read_json(folder/'partition.json');grid=read_json(work/'grid.json')
    transform=rasterio.Affine(*grid['transform'])
    database=folder/'models.sqlite'
    if sha256(database)!=fitted['sha256']:raise ValueError('Completed national profiles changed')
    for name,digest in partitioned['files'].items():
        if sha256(folder/name)!=digest:raise ValueError('National geometry index changed')
    ids=np.load(folder/'nodes.npy',mmap_mode='r');labels=np.load(folder/'labels.npy',mmap_mode='r')
    outdir=root/'output'/label/'cogs';outdir.mkdir(parents=True,exist_ok=True)
    recipe=digest_json(dict(fit=fitted['sha256'],code=sha256(Path(__file__)),
                           cog_code=sha256(Path(__file__).parents[1]/'cog.py')))
    reports=[];audits=[];joins=[];borders=[]
    with sqlite3.connect(f'file:{database}?mode=ro',uri=True) as db, \
         rasterio.Env(GDAL_CACHEMAX=48*1024**2,GDAL_NUM_THREADS='1'),rasterio.open(work/'classes.vrt') as mask:
        db.execute('PRAGMA cache_size=-32768')
        for number,record in enumerate(geometry['sheets']):
            name=record['sheet'];source=root/'sheets/coastal'/name;output=outdir/name
            source_sha=geometry['recipe']['sheets'][name]
            if sha256(source)!=source_sha:raise ValueError('Prepared sheet changed before river application')
            path=work/'geometry'/(Path(name).stem+'.npz')
            if sha256(path)!=record['sha256']:raise ValueError('Geometry cache changed before river application')
            with np.load(path) as data:pixels,owners=data['pixels'],data['owners']
            candidate,groups=project_owners(owners,ids,labels,db)
            with rasterio.open(source) as src:
                before=src.read(1);profile=src.profile.copy();tags=src.tags();t=src.transform
                x,y=(~transform)*(t.c,t.f);x,y=round(x),round(y)
                classes=mask.read(1,window=Window(x,y,src.width,src.height))
            row=pixels//grid['width']-y;col=pixels%grid['width']-x
            if np.any((row<0)|(row>=before.shape[0])|(col<0)|(col>=before.shape[1])):
                raise ValueError('River projection escaped its sheet')
            valid=(before!=-255)&np.isfinite(before)
            treated=np.isfinite(candidate)&valid[row,col]
            if np.any(candidate[treated]==-255):raise ValueError('River height collides with nodata')
            expected=before.copy();expected[row[treated],col[treated]]=candidate[treated]
            report_path=work/(Path(name).stem+'-cog.json')
            if output.exists():
                report=read_json(report_path)
                if report['recipe_sha256']!=recipe or report['sha256']!=sha256(output):
                    raise ValueError('Existing national COG changed or has a different recipe')
            else:
                temporary=root/'tmp'/('river-'+name)
                try:
                    with rasterio.open(temporary,'w',**profile) as dst:
                        dst.write(expected,1);dst.update_tags(**tags)
                        dst.update_tags(RIVER_TREATMENT='Shared connected CartoDEM profiles; unsupported components unchanged',
                                        RIVER_TREATED_PIXELS=int(treated.sum()),RIVER_PROFILE_SHA256=fitted['sha256'],
                                        RIVER_MASK=geometry['recipe']['river_source'],
                                        NONRIVER_ELEVATIONS='Unchanged',VERTICAL_DATUM='EGM2008',
                                        QUALITY_STATUS='Candidate build; native gradient diagnostics require review')
                    report=export_cog(temporary,output)
                finally:temporary.unlink(missing_ok=True)
                delta=candidate[treated]-before[row[treated],col[treated]]
                report.update(recipe_sha256=recipe,source_sha256=source_sha,river_pixels=len(pixels),
                              valid_river_pixels=int(valid[row,col].sum()),treated_river_pixels=int(treated.sum()),
                              nodata_pixels=int((~valid).sum()),
                              p95_absolute_river_change_m=float(np.percentile(np.abs(delta),95)) if len(delta) else None,
                              maximum_absolute_river_change_m=float(np.abs(delta).max()) if len(delta) else None)
                atomic_json(report_path,report)
            with rasterio.open(output) as dst:
                after=dst.read(1)
                if dst.nodata!=-255 or dst.transform!=t or dst.shape!=before.shape:
                    raise ValueError('Final COG geometry/nodata differ')
                if not np.array_equal(after,expected,equal_nan=True):raise ValueError('Final COG differs from shared profile projection')
                if not np.array_equal(dst.read_masks(1)>0,valid):raise ValueError('Final COG changes native validity')
            river=valid&(classes==3)
            if np.any(expected[~river]!=before[~river]):raise ValueError('Nonriver or nodata pixels changed')
            if np.any(after[classes==1]!=-255):raise ValueError('Ocean-mask pixels are not nodata')
            original=before[row,col];actual=after[row,col]
            check=audit_native_joins(pixels,original,actual,treated,groups,grid['width'])
            if check['boundary_pairs']:raise ValueError('Partial river treatment interface inside sheet')
            joins.append(dict(sheet=name,**check))
            border=(row==0)|(col==0)|(row==before.shape[0]-1)|(col==before.shape[1]-1)
            borders.append(tuple(a[border] for a in (pixels,original,actual,treated,groups))+
                           (np.full(int(border.sum()),number,np.int32),))
            item=dict(sheet=name,output_sha256=sha256(output),valid_river_pixels=int(river.sum()),
                      changed_river_pixels=int(np.sum(river&(before!=after))),native_validity_unchanged=True,
                      nonriver_and_nodata_unchanged=True,ocean_is_nodata=True,
                      adjacent_water_jumps=gradient_report(before,after,river,t))
            atomic_json(work/(Path(name).stem+'-audit.json'),item)
            reports.append(report);audits.append(item)
            emit('river_cog',complete=number+1,total=len(geometry['sheets']),sheet=name,
                 treated=report['treated_river_pixels'],valid_river=report['valid_river_pixels'])
    if borders:
        p,b,a,t,g,s=(np.concatenate([r[i] for r in borders]) for i in range(6))
        check=audit_native_joins(p,b,a,t,g,grid['width'],s)
        if check['boundary_pairs']:raise ValueError('Partial river treatment interface across sheet boundary')
        joins.append(dict(sheet='Across sheet borders',**check))
    audit=dict(label=label,partitioned=True,auditor_sha256=sha256(Path(__file__)),sheets=audits,
               downstream_profile_violations=0,native_treatment_joins=joins,
               fit_rejections=fitted['fit_rejections'],requires_gradient_review=True,
               maximum_added_water_jump_m=max((a['maximum_added_jump_m'] for s in audits for a in s['adjacent_water_jumps']),default=0.))
    atomic_json(work/'audit.json',audit)
    summary=dict(recipe_sha256=recipe,sheets=reports,nodata=-255,
                 river_source=geometry['recipe']['river_source'],
                 treated_river_pixels=sum(r['treated_river_pixels'] for r in reports),
                 valid_river_pixels=sum(r['valid_river_pixels'] for r in reports),
                 nodata_pixels=sum(r['nodata_pixels'] for r in reports),
                 quality_status='Candidate build; gradient review precedes publication')
    atomic_json(work/'cogs.json',summary)
    return summary


def main():
    import argparse
    import fcntl
    from build import setup
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['routing','fit','apply'])
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--label',default='final')
    args=parser.parse_args()
    if Path(args.label).name!=args.label:parser.error('Use a plain output label')
    root=setup(args.workdir)
    with (root/'logs/river-build.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if args.stage=='routing':
            from .regional import route_regions
            route_regions(root,args.label)
        elif args.stage=='fit':fit(root,args.label)
        else:apply(root,args.label)


if __name__=='__main__':main()
