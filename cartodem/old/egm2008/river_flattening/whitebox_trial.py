"""Measure Whitebox routing and watershed delineation on an analysis DEM crop.

This is an isolated comparison, not a production backend switch. The conditioned
DEM is disposable routing input and never becomes a published elevation sheet.
"""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np
import rasterio
from rasterio.windows import Window

from build import atomic_json, read_json, sha256
from .routing import analysis_window
from .whitebox import TO_DIRECTION, uncertain_upstream, measured, grid_matches, compatible_tiff, installation, run_flow


def prepare_context(root, folder, geometry_label, sheet, halo=1024):
    """Prepare the same overlapping DEM context without invoking any GRASS tool."""
    work = root/'rivers'/geometry_label
    with rasterio.Env(GDAL_CACHEMAX=48*1024**2, GDAL_NUM_THREADS='1', GDAL_MAX_DATASET_POOL_SIZE=24), \
            rasterio.open(work/'coastal.vrt') as src, rasterio.open(work/'classes.vrt') as classes, \
            rasterio.open(root/'sheets/coastal'/sheet) as core:
        x,y=(~src.transform)*(core.transform.c,core.transform.f)
        x,y=round(x),round(y)
        x0,y0=max(0,x-halo),max(0,y-halo)
        x1,y1=min(src.width,x+core.width+halo),min(src.height,y+core.height+halo)
        width,height=x1-x0,y1-y0
        if width*height>40_000_000:
            raise ValueError('Trial context exceeds 40 million cells')
        profile=dict(driver='GTiff',width=width,height=height,crs=src.crs,count=1,
                     transform=src.window_transform(Window(x0,y0,width,height)),
                     tiled=True,blockxsize=256,blockysize=256,compress='DEFLATE')
        with rasterio.open(folder/'analysis.tif','w',**profile,dtype='float32',nodata=-255,predictor=3) as dst, \
                rasterio.open(folder/'evidence.tif','w',**profile,dtype='uint8') as ev:
            for yy in range(0,height,512):
                for xx in range(0,width,512):
                    w,h=min(512,width-xx),min(512,height-yy)
                    a,b=max(0,xx-72),max(0,yy-72)
                    c,d=min(width,xx+w+72),min(height,yy+h+72)
                    window=Window(x0+a,y0+b,c-a,d-b)
                    values,valid,_=analysis_window(src.read(1,window=window),classes.read(1,window=window))
                    crop=np.s_[yy-b:yy-b+h,xx-a:xx-a+w]
                    dst.write(values[crop],1,window=Window(xx,yy,w,h))
                    ev.write(valid[crop].astype('uint8'),1,window=Window(xx,yy,w,h))


def sample(path, rows, cols):
    y, x = int(rows.min()), int(cols.min())
    window = Window(x, y, int(cols.max()-x+1), int(rows.max()-y+1))
    with rasterio.open(path) as ds:
        return ds.read(1, window=window)[rows-y, cols-x]


def compare(root, folder, geometry_label, sheet):
    with rasterio.open(folder/'analysis.tif') as src, rasterio.open(folder/'pointer.tif') as pntr:
        if not grid_matches(src, pntr):
            raise ValueError('Whitebox pointer grid changed')
        valid = src.read_masks(1) > 0
        values = pntr.read(1)
        if np.any(valid & ~np.isin(values, [0, 1, 2, 4, 8, 16, 32, 64, 128])):
            raise ValueError('Unexpected Whitebox pointer value')
        if not np.array_equal(pntr.read_masks(1)>0, valid):
            raise ValueError('Whitebox pointer nodata footprint changed')
        pointer = np.where(valid, values, 0).astype('uint8')
        del values
        profile = src.profile.copy()
        transform = src.transform
        for name in ('conditioned.tif', 'accumulation.tif'):
            with rasterio.open(folder/name) as other:
                if not grid_matches(src, other):
                    raise ValueError(f'{name}: grid changed')
        uncertainty = uncertain_upstream(pointer, valid)
        profile.update(dtype='uint8', nodata=None, predictor=1)
        with rasterio.open(folder/'uncertain.tif', 'w', **profile) as dst:
            dst.write(uncertainty.astype('uint8'), 1)

    work = root/'rivers'/geometry_label
    grid = read_json(work/'grid.json')
    native = rasterio.Affine(*grid['transform'])
    x, y = (~native)*(transform.c, transform.f)
    x, y = round(x), round(y)
    with np.load(work/'geometry'/(Path(sheet).stem+'.npz')) as geometry:
        pixels = geometry['pixels']
    rows, cols = pixels//grid['width']-y, pixels%grid['width']-x
    if (np.any(rows<0) or np.any(cols<0) or np.any(rows>=valid.shape[0])
            or np.any(cols>=valid.shape[1])):
        raise ValueError('Saved context does not cover this sheet geometry')
    evidence = sample(folder/'evidence.tif', rows, cols)>0
    evidence &= (rows>=8)&(cols>=8)&(rows<valid.shape[0]-8)&(cols<valid.shape[1]-8)
    wdir = TO_DIRECTION[pointer[rows, cols]]
    warea = sample(folder/'accumulation.tif', rows, cols)
    wok = evidence & ~uncertainty[rows, cols] & (wdir>0) & (warea>0) & np.isfinite(warea)
    result = dict(river_pixels=len(pixels), original_evidence_pixels=int(evidence.sum()),
                  whitebox_accepted_pixels=int(wok.sum()),
                  valid_analysis_pixels=int(valid.sum()), uncertain_pixels=int(uncertainty.sum()))
    if (folder/'grass-drainage.tif').exists():
        gdir = sample(folder/'grass-drainage.tif', rows, cols)
        garea = sample(folder/'grass-accumulation.tif', rows, cols)
        gok = evidence & (gdir>0) & (gdir<=8) & (garea>0) & np.isfinite(garea)
        common = wok & gok
        difference = np.abs(wdir[common].astype(int)-gdir[common].astype(int))
        turns = np.minimum(difference, 8-difference)
        result['grass_comparison'] = dict(accepted_pixels=int(gok.sum()), common_pixels=int(common.sum()),
            direction_difference_degrees={str(i*45):int((turns==i).sum()) for i in range(5)},
            whitebox_only_pixels=int((wok&~gok).sum()), grass_only_pixels=int((gok&~wok).sum()))
    changed = 0
    delta_min = delta_max = 0.
    with rasterio.open(folder/'analysis.tif') as src, rasterio.open(folder/'conditioned.tif') as conditioned:
        for _, window in src.block_windows(1):
            a = src.read(1, window=window)
            b = conditioned.read(1, window=window)
            good = src.read_masks(1, window=window)>0
            if not np.array_equal(good, conditioned.read_masks(1, window=window)>0):
                raise ValueError('Routing conditioning changed nodata footprint')
            delta = b[good]-a[good]
            changed += int(np.count_nonzero(delta))
            if len(delta):
                delta_min = min(delta_min, float(delta.min()))
                delta_max = max(delta_max, float(delta.max()))
    result['disposable_conditioning'] = dict(changed_pixels=changed, minimum_delta_m=delta_min,
                                           maximum_delta_m=delta_max, published=False)
    if not wok.any():
        raise ValueError('No reliable WRIS river sample for the watershed check')
    index = np.flatnonzero(wok)[np.argmax(warea[wok])]
    r, c = int(rows[index]), int(cols[index])
    lon, lat = transform*(c+.5, r+.5)
    result['pour_point'] = dict(row=r, col=c, latitude=lat, longitude=lon,
                               expected_upstream_cells=float(warea[index]))
    profile.update(dtype='int32', nodata=0, predictor=2)
    with rasterio.open(folder/'pour-point.tif', 'w', **profile) as dst:
        dst.write(np.ones((1,1), np.int32), 1, window=Window(c,r,1,1))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--context', type=Path,
                        help='Optional existing analysis.tif/evidence.tif directory; otherwise prepare from corrected sheets')
    parser.add_argument('--grass-reference', type=Path,
                        help='Optional existing drainage.tif/accumulation.tif directory for comparison only')
    parser.add_argument('--sheet', required=True)
    parser.add_argument('--geometry-label', default='final')
    parser.add_argument('--label', default='whitebox-routing-v1')
    parser.add_argument('--max-rss-mib', type=int, default=1100)
    parser.add_argument('--workflow', choices=('staged','full'), default='staged',
                        help='Separate breaching/D8 processes, or the aspect-aligned combined workflow')
    parser.add_argument('--tool-installation',default='tools/whitebox-prebuilt',
                        help='Workspace-relative directory containing the tool build.json')
    args = parser.parse_args()
    if Path(args.label).name!=args.label or args.label in ('.','..'):
        parser.error('Use a plain label')
    root = args.workdir.expanduser().resolve()
    folder = root/'rivers'/args.label
    folder.mkdir(exist_ok=False)
    build,binary,tool_env=installation(root,args.tool_installation)
    report = dict(state='running', tool=build, sheet=args.sheet,
                  geometry_label=args.geometry_label, production_backend_changed=False,
                  published_heights_changed=False, workflow=args.workflow,
                  method=('BreachDepressions + D8Pointer + D8FlowAccumulation --pntr; cells' if args.workflow=='staged'
                          else 'FlowAccumulationFullWorkflow; cells; default pointer'),
                  code_sha256=sha256(Path(__file__)))
    try:
        hashes = {}
        inputs = []
        if args.context:
            inputs.extend([(args.context/'analysis.tif', 'analysis.tif'),
                           (args.context/'evidence.tif', 'evidence.tif')])
        else:
            prepare_context(root, folder, args.geometry_label, args.sheet)
            hashes.update((name,sha256(folder/name)) for name in ('analysis.tif','evidence.tif'))
        if args.grass_reference:
            inputs.extend((args.grass_reference/source, target) for source,target in [
                ('drainage.tif','grass-drainage.tif'), ('accumulation.tif','grass-accumulation.tif')])
        for source, target in inputs:
            shutil.copy2(source, folder/target)
            hashes[target] = sha256(folder/target)
        if args.grass_reference and (args.grass_reference/'analysis.tif').exists():
            report['same_analysis_as_reference'] = hashes['analysis.tif']==sha256(args.grass_reference/'analysis.tif')
            if not report['same_analysis_as_reference']:
                raise ValueError('Comparison requires the exact same analysis raster')
        report['inputs_sha256'] = hashes
        atomic_json(folder/'report.json', report)
        compatible_tiff(folder)
        hashes['analysis-whitebox.tif']=sha256(folder/'analysis-whitebox.tif')
        report['whitebox_input_encoding']=dict(predictor=1,pixel_values_verified_identical=True)
        report['routing_runs']={}
        for name,metrics in run_flow(binary,folder,tool_env,workflow=args.workflow,max_rss_mib=args.max_rss_mib):
            report['routing_runs'][name]=metrics
            atomic_json(folder/'report.json',report)
        with rasterio.Env(GDAL_CACHEMAX=48*1024**2, GDAL_NUM_THREADS='1'):
            report['comparison'] = compare(root, folder, args.geometry_label, args.sheet)
        atomic_json(folder/'report.json', report)
        report['watershed_run'] = measured([
            str(binary), '-r=Watershed', '-v', f'--d8_pntr={folder/"pointer.tif"}',
            f'--pour_pts={folder/"pour-point.tif"}', f'--output={folder/"watershed.tif"}',
        ], folder/'watershed.log', args.max_rss_mib,tool_env)
        with rasterio.open(folder/'watershed.tif') as watershed, rasterio.open(folder/'analysis.tif') as src:
            if not grid_matches(src, watershed):
                raise ValueError('Watershed grid changed')
            cells = sum(int(np.count_nonzero(watershed.read(1,window=w)==1)) for _,w in watershed.block_windows(1))
        expected = report['comparison']['pour_point']['expected_upstream_cells']
        report['watershed_cells'] = cells
        report['watershed_matches_accumulation'] = bool(np.float32(cells)==np.float32(expected))
        if not report['watershed_matches_accumulation']:
            raise ValueError('Watershed area disagrees with flow accumulation')
        for name, digest in hashes.items():
            if sha256(folder/name)!=digest:
                raise ValueError(f'Comparison input changed: {name}')
        report['outputs_sha256'] = {name:sha256(folder/name) for name in (
            'conditioned.tif','pointer.tif','accumulation.tif','uncertain.tif','pour-point.tif','watershed.tif')}
        report['state'] = 'complete'
    except BaseException as exc:
        report.update(state='failed', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        atomic_json(folder/'report.json', report)
        print(json.dumps(report), flush=True)


if __name__=='__main__':
    main()
