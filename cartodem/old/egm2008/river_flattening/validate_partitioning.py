"""Compare partitioned and whole-network fits with identical routing evidence.

Run only on bounded validation domains; the whole-network comparison deliberately
holds all pilot observations in memory. Published COGs are never modified.
"""
import argparse
from pathlib import Path
import sqlite3

import numpy as np
import rasterio

from build import atomic_json,read_json,sha256
from .boundaries import contact_files,native_values
from .geometry import compress_network
from .joint import fit_joint
from .national import lookup,project_owners


def validate(root,label):
    root=Path(root);work=root/'rivers'/label;records=read_json(work/'geometry.json')
    grid=read_json(work/'grid.json')
    if grid['width']*grid['height']>40_000_000:
        raise ValueError('Validate equivalence only on a bounded pilot')
    collected={k:[] for k in ('nodes','bank_owners','bank_heights','ocean_owners')}
    for record in records['sheets']:
        path=work/'geometry'/(Path(record['sheet']).stem+'.npz')
        if sha256(path)!=record['sha256']:raise ValueError('Geometry cache changed')
        with np.load(path) as data:
            for k in collected:collected[k].append(data[k])
    native={k:np.concatenate(v) for k,v in collected.items()};del collected
    native['nodes']=np.unique(native['nodes']);nodes=native['nodes']
    routing=dict(vector=np.zeros((len(nodes),2)),weight=np.zeros(len(nodes)),count=np.zeros(len(nodes),np.int64))
    for record in read_json(work/'routing/report.json')['regions']:
        path=work/'routing'/record['evidence']
        if sha256(path)!=record['sha256']:raise ValueError('Routing evidence changed')
        with np.load(path) as data:
            index=lookup(nodes,data['nodes'])
            for axis in range(2):np.add.at(routing['vector'][:,axis],index,data['vector'][:,axis])
            for k in ('weight','count'):np.add.at(routing[k],index,data[k])
    t=rasterio.Affine(*grid['transform']);lat=t.f+t.e*(nodes//grid['width']+.5)
    network=compress_network(nodes,grid['width'],t.a*111195.0802*np.cos(np.deg2rad(lat)),abs(t.e)*111195.0802)
    network,_=fit_joint(native,network,routing=routing,width=grid['width'],contacts=contact_files(root,work))
    expected,_=native_values(network)
    folder=work/'components';ids=np.load(folder/'nodes.npy');labels=np.load(folder/'labels.npy')
    with sqlite3.connect(f'file:{folder/"models.sqlite"}?mode=ro',uri=True) as db:
        actual,_=project_owners(nodes,ids,labels,db)
    if not np.array_equal(np.isfinite(actual),np.isfinite(expected)):
        raise ValueError('Partitioning changes the set of treated skeleton nodes')
    valid=np.isfinite(expected);delta=np.abs(actual[valid]-expected[valid])
    maximum=float(delta.max(initial=0))
    if maximum>.001:raise ValueError(f'Partitioning changes fitted heights by {maximum} m')
    report=dict(native_nodes=len(nodes),treated_skeleton_nodes=int(valid.sum()),
                identical_treatment=True,maximum_height_difference_m=maximum,tolerance_m=.001,
                geometry_sha256=sha256(work/'geometry.json'),routing_sha256=sha256(work/'routing/report.json'),
                fit_sha256=sha256(work/'fit.json'),validation_code_sha256=sha256(Path(__file__)))
    atomic_json(work/'partition-equivalence.json',report)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,required=True);parser.add_argument('--label',required=True)
    args=parser.parse_args();print(validate(args.workdir,args.label),flush=True)


if __name__=='__main__':main()
