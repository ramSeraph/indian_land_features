"""Compare production owner evidence with independently generated trial rasters."""
import argparse
from pathlib import Path

import numpy as np
import rasterio

from build import atomic_json, read_json, sha256
from .regional import accumulate
from .whitebox import TO_DIRECTION
from .whitebox_trial import sample


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--label',default='final')
    parser.add_argument('--trial',default='whitebox-prebuilt-v1')
    parser.add_argument('--sheet',default='cdnc43d.tif')
    args=parser.parse_args();root=args.workdir.expanduser().resolve()
    work=root/'rivers'/args.label;trial=root/'rivers'/args.trial
    report=read_json(trial/'report.json')
    if report['state']!='complete' or report['geometry_label']!=args.label or report['sheet']!=args.sheet:
        raise ValueError('A completed trial using the same geometry is required')
    for name,digest in {**report['inputs_sha256'],**report['outputs_sha256']}.items():
        if sha256(trial/name)!=digest:raise ValueError('Trial raster changed: '+name)
    grid=read_json(work/'grid.json')
    with rasterio.open(trial/'analysis.tif') as ds:
        x,y=(~rasterio.Affine(*grid['transform']))*(ds.transform.c,ds.transform.f)
        x,y=round(x),round(y);height,width=ds.shape
    stem=Path(args.sheet).stem
    with np.load(work/'geometry'/(stem+'.npz')) as data:
        pixels,owners=data['pixels'],data['owners']
    rows,cols=pixels//grid['width']-y,pixels%grid['width']-x
    if np.any((rows<0)|(cols<0)|(rows>=height)|(cols>=width)):
        raise ValueError('Trial crop does not cover the owning sheet')
    with rasterio.Env(GDAL_CACHEMAX=32*1024**2):
        valid=sample(trial/'evidence.tif',rows,cols)>0
        valid &= sample(trial/'uncertain.tif',rows,cols)==0
        valid &= (rows>=8)&(cols>=8)&(rows<height-8)&(cols<width-8)
        values=sample(trial/'pointer.tif',rows,cols)
        direction=np.zeros(values.shape,np.uint8)
        direction[valid]=TO_DIRECTION[values[valid].astype('uint8')]
        area=sample(trial/'accumulation.tif',rows,cols)
    expected=accumulate(owners,direction,area,valid)
    target=work/'routing'/(stem+'-owners.npz')
    checkpoint=read_json(target.with_suffix('.json'))
    if sha256(target)!=checkpoint['sha256']:raise ValueError('Production evidence checksum changed')
    if checkpoint['temporary_raster_sha256']['accumulation.tif']!=report['outputs_sha256']['accumulation.tif']:
        raise ValueError('Production accumulation raster differs from the trial')
    with np.load(target) as actual:
        for name,values in expected.items():
            np.testing.assert_array_equal(actual[name],values,err_msg=name)
    result=dict(state='complete',sheet=args.sheet,label=args.label,trial=args.trial,
                sampled_pixels=int(expected['count'].sum()),owners=len(expected['nodes']),
                evidence_arrays_exactly_equal=True,accumulation_raster_sha256_equal=True,
                evidence_sha256=sha256(target),trial_report_sha256=sha256(trial/'report.json'),
                code_sha256=sha256(Path(__file__)),production_context=checkpoint)
    output=work/'routing'/(stem+'-trial-equivalence.json')
    atomic_json(output,result)
    print(f'Exact trial/production evidence match: {result["sampled_pixels"]} samples; {result["owners"]} owners',flush=True)


if __name__=='__main__':main()
