"""Bounded river pilot: extract geometry, fit supported networks, audit changes.

This is a study runner, not the nationwide production stage. It makes direction
candidates explicit and leaves unsupported components unchanged. All outputs
are new files; neither input raster is edited.
"""
import argparse
import hashlib
import importlib.metadata
import json
import math
import shutil
from pathlib import Path

import numpy as np
import rasterio
from scipy import sparse
from scipy.sparse.csgraph import connected_components, dijkstra
from rasterio.windows import Window, from_bounds

from .geometry import extract_window, compress_network, orient_component
from .profiles import fit_network


def assign_stations(network, owners):
    """Map native skeleton owners to the nearest longitudinal station."""
    ids = network['ids']
    indices = np.searchsorted(ids, owners)
    safe = np.minimum(indices,len(ids)-1)
    good = (indices<len(ids)) & (ids[safe]==owners)
    stations = np.where(network['fraction'][safe]<=.5,network['left'][safe],network['right'][safe])
    good &= stations>=0
    return stations,good


def estimate_surface(height, classes, *, dx, dy, spacing=180., quantile=.2,
                     stiffness=2., minimum_outlet_drop=3., max_support_gap=600.):
    """Return pilot elevations, detailed component reports and fitted geometry."""
    native = extract_window(height,classes,dx=dx,dy=dy,grid_width=height.shape[1])
    if not len(native['nodes']):
        return height.copy(),dict(components=[],river_pixels=0,modified_pixels=0),None
    network = compress_network(native['nodes'],height.shape[1],dx,dy,spacing=spacing)
    network,reports=fit_geometry(native,network,spacing=spacing,quantile=quantile,stiffness=stiffness,
                                 minimum_outlet_drop=minimum_outlet_drop,max_support_gap=max_support_gap)
    return apply_surface(height,classes,native,network,reports)


def fit_geometry(native,network,*,spacing=180.,quantile=.2,stiffness=2.,minimum_outlet_drop=3.,max_support_gap=600.):
    """Fit assembled geometry once, independently of sheet boundaries."""
    n = len(network['station_pixels'])
    obs_station,obs_ok = assign_stations(network,native['bank_owners'])
    obs_station,obs_height = obs_station[obs_ok],native['bank_heights'][obs_ok]
    coast_station,coast_ok = assign_stations(network,native['ocean_owners'])
    coast_station = np.unique(coast_station[coast_ok])
    edges = {}
    for reach,x in zip(network['reaches'],network['distances']):
        for u,v,d in zip(reach[:-1],reach[1:],np.diff(x)):
            key=tuple(sorted((u,v)))
            edges[key]=min(edges.get(key,np.inf),float(d))
    edge_pairs=np.asarray(list(edges),dtype=int).reshape(-1,2)
    if not len(edge_pairs):
        network.update(fitted=np.full(n,np.nan,dtype=np.float32),supported_edges=np.empty((0,2),dtype=int),
                       bank_reference=np.full(n,np.nan),bank_count=np.bincount(obs_station,minlength=n),
                       component=np.arange(n),missing_bank_owners=int((~obs_ok).sum()))
        return network,[]
    lengths=np.asarray(list(edges.values()))
    graph=sparse.csr_matrix((np.tile(lengths,2),
                            (np.r_[edge_pairs[:,0],edge_pairs[:,1]],np.r_[edge_pairs[:,1],edge_pairs[:,0]])),shape=(n,n))
    nc,labels=connected_components(graph,directed=False)
    count=np.bincount(obs_station,minlength=n)
    reference=np.full(n,np.nan)
    order=np.argsort(obs_station,kind='stable')
    starts=np.r_[0,np.cumsum(count)]
    for i in np.flatnonzero(count>=4):
        reference[i]=np.quantile(obs_height[order[starts[i]:starts[i+1]]],quantile)
    fitted=np.full(n,np.nan,dtype=np.float32)
    reports=[]
    supported_edges=[]
    # Assemble component slices once; avoid scanning all observations/reaches for
    # every small disconnected river in the nationwide inventory.
    member_order=np.argsort(labels,kind='stable')
    member_starts=np.r_[0,np.cumsum(np.bincount(labels,minlength=nc))]
    obs_labels=labels[obs_station]
    obs_order=np.argsort(obs_labels,kind='stable')
    obs_starts=np.r_[0,np.cumsum(np.bincount(obs_labels,minlength=nc))]
    reach_groups=[[] for _ in range(nc)]
    for reach,x in zip(network['reaches'],network['distances']):
        reach_groups[labels[reach[0]]].append((reach,x))
    for component in range(nc):
        members=member_order[member_starts[component]:member_starts[component+1]]
        report=dict(component=int(component),stations=len(members),status='unchanged')
        reports.append(report)
        local_graph=graph[members][:,members]
        terminals=np.flatnonzero(np.diff(local_graph.indptr)==1)
        if len(members)<3 or len(terminals)<2:
            report['reason']='Too few stations or no unambiguous terminal pair'; continue
        observed=np.flatnonzero(np.isfinite(reference[members]))
        report['supported_stations']=len(observed)
        if len(observed)<3 or len(observed)/len(members)<.6:
            report['reason']='Insufficient CartoDEM bank support'; continue
        distance,_,nearest=dijkstra(local_graph,directed=False,indices=observed,min_only=True,return_predecessors=True)
        report['maximum_support_gap_m']=float(distance.max())
        if distance.max()>max_support_gap:
            report['reason']='A section is too far from supported bank observations'; continue
        local_reference=reference[members].copy()
        filled_reference=reference[members[nearest]]
        coastal=np.flatnonzero(np.isin(members,coast_station))
        coast_terminals=np.array([],dtype=int)
        if len(coastal):
            coast_dist=dijkstra(local_graph,directed=False,indices=coastal,min_only=True)
            coast_terminals=terminals[coast_dist[terminals]<=2*spacing]
        terminal_levels=filled_reference[terminals]
        report['terminal_candidates']=[dict(station=int(members[t]),height_m=float(filled_reference[t]),
                                                  near_ocean=bool(t in coast_terminals)) for t in terminals]
        if len(coast_terminals)==1:
            root=int(coast_terminals[0]); evidence='Single terminal adjoining the GLO-30 ocean class'
        elif len(coast_terminals)>1:
            report['reason']='Multiple ocean outlets require explicit direction review'; continue
        else:
            ranked=np.argsort(terminal_levels)
            separation=float(terminal_levels[ranked[1]]-terminal_levels[ranked[0]])
            report['outlet_height_separation_m']=separation
            if separation<minimum_outlet_drop:
                report['reason']='Terminal bank levels do not establish a distinct low outlet'; continue
            root=int(terminals[ranked[0]]); evidence='Distinct lowest terminal bank quantile (direction hypothesis)'
        report['outlet_station']=int(members[root])
        report['direction_evidence']=evidence
        remap=np.full(n,-1,dtype=int); remap[members]=np.arange(len(members))
        reaches=[]; distances=[]
        for reach,x in reach_groups[component]:
            reaches.append(remap[reach].tolist()); distances.append(x)
        try:
            directed,rr,dd,potential=orient_component(len(members),reaches,distances,root)
            # Direction must agree with supported CartoDEM bank levels on long
            # reaches. A harmonic potential alone is not evidence of real flow.
            disagreements=[]
            for r,x in zip(rr,dd):
                upstream,downstream=filled_reference[r[0]],filled_reference[r[-1]]
                if downstream-upstream>5:
                    disagreements.append(dict(upstream=int(members[r[0]]),downstream=int(members[r[-1]]),
                                              uphill_bank_difference_m=float(downstream-upstream)))
            if disagreements:
                report['direction_conflicts']=disagreements
                if len(coast_terminals)!=1:
                    raise ValueError('Direction hypothesis conflicts with bank elevations')
                report['direction_conflict_policy']='Ocean outlet supplies direction; bank disagreements remain fit diagnostics'
            use=obs_order[obs_starts[component]:obs_starts[component+1]]
            upper=np.where(np.isfinite(local_reference),local_reference,np.inf)
            values,diagnostic=fit_network(len(members),directed,remap[obs_station[use]],obs_height[use],
                                          reaches=rr,distances=dd,quantile=quantile,
                                          stiffness=stiffness,upper_bounds=upper)
            difference=np.abs(values[observed]-local_reference[observed])
            report['p95_station_adjustment_m']=float(np.percentile(difference,95))
            report['maximum_station_adjustment_m']=float(difference.max())
            if np.percentile(difference,95)>15 or difference.max()>40:
                raise ValueError('Profile requires excessive adjustment of supported bank levels')
            fitted[members]=values
            supported_edges.extend(members[directed].tolist())
            report.update(status='fitted',fit=diagnostic)
        except ValueError as exc:
            report['reason']=str(exc)
    network.update(fitted=fitted,supported_edges=np.asarray(supported_edges),
                   bank_reference=reference,bank_count=count,component=labels,
                   missing_bank_owners=int((~obs_ok).sum()))
    return network,reports


def apply_surface(height,classes,native,network,reports):
    fitted=network['fitted']
    # Project along the in-water owner assignment. This never writes a land,
    # lake, ocean or original nodata pixel.
    owner_index=np.searchsorted(network['ids'],native['owners'])
    left,right=network['left'][owner_index],network['right'][owner_index]
    fraction=network['fraction'][owner_index]
    good=(left>=0)&(right>=0)
    candidate=np.full(len(left),np.nan,dtype=np.float32)
    candidate[good]=(1-fraction[good])*fitted[left[good]]+fraction[good]*fitted[right[good]]
    original=height.ravel()[native['pixels']]
    valid=np.isfinite(original)&(original!=-255)
    apply=valid&np.isfinite(candidate)
    if np.any(candidate[apply]==-255):
        raise ValueError('River surface collides with the reserved nodata sentinel')
    output=height.copy()
    output.ravel()[native['pixels'][apply]]=candidate[apply]
    outside=(classes!=3)|~np.isfinite(height)|(height==-255)
    if not np.array_equal(output[outside],height[outside],equal_nan=True):
        raise AssertionError('River treatment changed land, lake, ocean or source nodata')
    delta=candidate[apply]-original[apply]
    report=dict(components=reports,river_pixels=len(native['pixels']),valid_river_pixels=int(valid.sum()),
                treated_river_pixels=int(apply.sum()),modified_pixels=int(np.count_nonzero(delta)),
                unsupported_skeleton_pixels=network['unsupported_skeleton_pixels'],
                missing_bank_owners=network['missing_bank_owners'],
                p95_absolute_change_m=float(np.percentile(np.abs(delta),95)) if len(delta) else None,
                max_absolute_change_m=float(np.abs(delta).max()) if len(delta) else None,
                land_lake_and_nodata_unchanged=True)
    return output,report,network


def run(dem_path,mask_path,bounds,outdir,*,spacing=180.,quantile=.2,stiffness=2.,max_support_gap=600.):
    outdir=Path(outdir)
    outdir.mkdir(parents=True,exist_ok=True)
    if (outdir/'report.json').exists():
        raise ValueError('Pilot output already exists; choose a new output directory')
    code_dir=outdir/'code'
    code_dir.mkdir(exist_ok=True)
    code_hashes={}
    for name in ('__init__.py','profiles.py','geometry.py','study.py'):
        source=Path(__file__).with_name(name)
        shutil.copyfile(source,code_dir/name)
        code_hashes[name]=hashlib.sha256(source.read_bytes()).hexdigest()
    with rasterio.Env(GDAL_CACHEMAX=64*1024**2,GDAL_NUM_THREADS='1'):
        with rasterio.open(dem_path) as dem,rasterio.open(mask_path) as mask:
            if dem.crs!=rasterio.crs.CRS.from_epsg(4326) or dem.crs!=mask.crs:
                raise ValueError('Expected EPSG:4326 elevation and classification rasters')
            # The pilot requires an already-aligned classification VRT/raster.
            if dem.shape!=mask.shape or not np.allclose(tuple(dem.transform),tuple(mask.transform),atol=1e-10,rtol=0):
                raise ValueError('Height and class rasters must share the native grid')
            floating=from_bounds(*bounds,transform=dem.transform)
            c0=max(0,math.floor(floating.col_off));r0=max(0,math.floor(floating.row_off))
            c1=min(dem.width,math.ceil(floating.col_off+floating.width));r1=min(dem.height,math.ceil(floating.row_off+floating.height))
            if (r1-r0)*(c1-c0)>20000000 or c1<=c0 or r1<=r0:
                raise ValueError('Pilot window is empty or exceeds the bounded-memory study limit')
            window=Window(c0,r0,c1-c0,r1-r0)
            h=dem.read(1,window=window).astype(np.float32)
            valid=(dem.read_masks(1,window=window)>0)&np.isfinite(h)
            if dem.nodata!=-255 and np.any(valid&(h==-255)):
                raise ValueError('A valid corrected elevation collides with the reserved -255 output sentinel')
            classes=mask.read(1,window=window)
            if np.any(classes>3):
                raise ValueError('Classification coverage is missing')
            valid &= classes!=1
            h[~valid]=-255
            transform=dem.window_transform(window)
            profile=dict(driver='GTiff',width=h.shape[1],height=h.shape[0],count=1,dtype='float32',
                         crs='EPSG:4326',transform=transform,nodata=-255.,tiled=True,
                         blockxsize=256,blockysize=256,compress='DEFLATE',predictor=3)
    latitude=transform.f+transform.e*h.shape[0]/2
    dy=abs(transform.e)*111195.0802
    dx=transform.a*111195.0802*math.cos(math.radians(latitude))
    output,report,network=estimate_surface(h,classes,dx=dx,dy=dy,spacing=spacing,quantile=quantile,stiffness=stiffness,
                                          max_support_gap=max_support_gap)
    for name,array in [('coastal_baseline',h),('river_pilot',output)]:
        with rasterio.open(outdir/f'{name}.tif','w',**profile) as dst:
            dst.write(array,1)
            dst.update_tags(VERTICAL_DATUM='EGM2008',STATUS='river method pilot; not production')
    with rasterio.open(outdir/'classes.tif','w',**{**profile,'dtype':'uint8','nodata':255,'predictor':1}) as dst:
        dst.write(classes,1)
    report.update(status='pilot',input_dem=Path(dem_path).name,input_mask=Path(mask_path).name,
                  code_sha256=code_hashes,
                  software={name:importlib.metadata.version(name) for name in ('numpy','scipy','scikit-image','clarabel','rasterio')},
                  input_window_sha256=dict(heights=hashlib.sha256(h.tobytes()).hexdigest(),
                                            classes=hashlib.sha256(classes.tobytes()).hexdigest()),
                  requested_bounds=list(bounds),transform=list(transform)[:6],shape=list(h.shape),
                  parameters=dict(station_spacing_m=spacing,quantile=quantile,stiffness=stiffness,
                                  bank_width_pixels=2,minimum_outlet_drop_m=3,max_support_gap_m=max_support_gap))
    features=[]
    if network is not None and 'fitted' in network:
        pixels=network['station_pixels']
        for i,pixel in enumerate(pixels):
            x,y=transform*(int(pixel%h.shape[1])+.5,int(pixel//h.shape[1])+.5)
            value=float(network['fitted'][i]);reference=float(network['bank_reference'][i])
            features.append(dict(type='Feature',geometry=dict(type='Point',coordinates=[x,y]),
                                 properties=dict(station=i,component=int(network['component'][i]),
                                                 fitted_m=value if math.isfinite(value) else None,
                                                 bank_quantile_m=reference if math.isfinite(reference) else None,
                                                 bank_samples=int(network['bank_count'][i]))))
    (outdir/'stations.geojson').write_text(json.dumps(dict(type='FeatureCollection',features=features),indent=2,allow_nan=False)+'\n')
    (outdir/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='components'},allow_nan=False),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dem',type=Path,required=True)
    p.add_argument('--mask',type=Path,required=True)
    p.add_argument('--bounds',nargs=4,type=float,required=True,metavar=('WEST','SOUTH','EAST','NORTH'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--station-spacing',type=float,default=180.)
    p.add_argument('--quantile',type=float,default=.2)
    p.add_argument('--stiffness',type=float,default=2.)
    p.add_argument('--max-support-gap',type=float,default=600.)
    a=p.parse_args()
    run(a.dem,a.mask,a.bounds,a.output,spacing=a.station_spacing,quantile=a.quantile,stiffness=a.stiffness,
        max_support_gap=a.max_support_gap)


if __name__=='__main__':
    main()
