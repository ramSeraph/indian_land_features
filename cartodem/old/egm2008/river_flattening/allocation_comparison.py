"""Compare allocation backends through their native-grid candidate COG outputs."""
import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio

from build import atomic_json,read_json,sha256
from .whitebox import grid_matches


def compare(root,label,previous):
    work=root/'rivers'/label;old=root/'rivers'/previous
    record=read_json(work/'geometry.json');rows=[]
    inputs={s['name']:s for s in read_json(root/'reports/inputs.json')['sheets']}
    for sheet in record['sheets']:
        name=sheet['sheet'];stem=Path(name).stem
        with np.load(work/'geometry'/(stem+'.npz')) as a,np.load(old/'geometry'/(stem+'.npz')) as b:
            for key in ('nodes','pixels','bank_heights','pruned_skeleton_pixels'):
                np.testing.assert_array_equal(a[key],b[key])
            owner_changes=int(np.sum(a['owners']!=b['owners']))
            bank_changes=int(np.sum(a['bank_owners']!=b['bank_owners']))
        before=root/'output'/previous/'cogs'/name;after=root/'output'/label/'cogs'/name
        changed=0;maximum=0.;nodata=0
        with rasterio.Env(GDAL_CACHEMAX=32*1024**2),rasterio.open(before) as a, \
             rasterio.open(after) as b,rasterio.open(root/'sheets/coastal'/name) as base, \
             rasterio.open(root/'inputs/glo30-wbm'/inputs[name]['mask']) as mask:
            for other in (a,base,mask):
                if not grid_matches(b,other):raise ValueError('Comparison grids differ')
            if a.nodata!=-255 or b.nodata!=-255:raise ValueError('Nodata sentinel changed')
            if b.tags().get('RIVER_MASK')!='glo30':raise ValueError('Unexpected river mask tag')
            for _,window in b.block_windows(1):
                x=a.read(1,window=window);y=b.read(1,window=window)
                original=base.read(1,window=window);classes=mask.read(1,window=window)
                invalid=original==-255
                np.testing.assert_array_equal(y==-255,invalid)
                np.testing.assert_array_equal(x==-255,invalid)
                preserve=(classes!=3)|invalid
                np.testing.assert_array_equal(y[preserve],original[preserve])
                changed+=int(np.sum(x!=y));nodata+=int(invalid.sum())
                maximum=max(maximum,float(np.max(np.abs(y-x))))
        rows.append(dict(sheet=name,water_owner_changes=owner_changes,bank_owner_changes=bank_changes,
                         changed_cog_pixels=changed,maximum_height_difference_m=maximum,nodata_pixels=nodata,
                         nodata_and_nonriver_preserved=True,before_sha256=sha256(before),after_sha256=sha256(after)))
    status=read_json(root/'reports'/f'build-{label}.json')
    report=dict(label=label,previous=previous,sheets=rows,
                elapsed_seconds=(datetime.fromisoformat(status['finished'])-datetime.fromisoformat(status['started'])).total_seconds(),
                cogs={k:v for k,v in read_json(work/'cogs.json').items() if k!='sheets'},
                audit={k:v for k,v in read_json(work/'audit.json').items() if k!='sheets'},
                code_sha256=sha256(Path(__file__)),publication_requires_gradient_review=True)
    atomic_json(work/'allocation-comparison.json',report)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',required=True,type=Path)
    parser.add_argument('--label',default='pilot-glo30-allocation-v1')
    parser.add_argument('--previous',default='pilot-glo30-whitebox-v1')
    args=parser.parse_args()
    print(compare(args.workdir.resolve(),args.label,args.previous))


if __name__=='__main__':main()
