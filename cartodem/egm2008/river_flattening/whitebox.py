"""Whitebox direction evidence; all conditioning stays on disposable rasters."""
import os
from pathlib import Path
import subprocess
import time

import numpy as np
import rasterio
from scipy.ndimage import binary_dilation

from build import atomic_json, read_json, sha256
from .routing import OFFSETS

TO_DIRECTION = np.zeros(256, np.uint8)
TO_DIRECTION[[1, 2, 4, 8, 16, 32, 64, 128]] = [1, 8, 7, 6, 5, 4, 3, 2]


def uncertain_upstream(pointer, valid):
    """Trace possible missing upstream area along the actual Whitebox pointers.

    Seed valid cells beside NULL or a crop edge, then mark their downstream
    paths. This is conservative at ocean edges. The supplied pointers define
    the network, including when evaluating an alternative routing algorithm.
    """
    if pointer.shape != valid.shape or pointer.ndim != 2:
        raise ValueError('Pointer and validity grids must match')
    if np.any(valid & (pointer != 0) & (TO_DIRECTION[pointer] == 0)):
        raise ValueError('Unknown Whitebox pointer code')
    seeds = binary_dilation(~valid, structure=np.ones((3, 3), bool), border_value=1) & valid
    height, width = pointer.shape
    marked = np.zeros(pointer.shape, bool)
    flat_marked, flat_pointer, flat_valid = marked.ravel(), pointer.ravel(), valid.ravel()
    for seed in np.flatnonzero(seeds):
        cell = int(seed)
        while not flat_marked[cell]:
            flat_marked[cell] = True
            direction = int(TO_DIRECTION[flat_pointer[cell]])
            if not direction:
                break
            r, c = divmod(cell, width)
            dr, dc = OFFSETS[direction]
            r, c = r + int(dr), c + int(dc)
            if not (0 <= r < height and 0 <= c < width):
                break
            cell = r * width + c
            if not flat_valid[cell]:
                break
    return marked


def measured(command, log_path, max_rss_mib, env=None):
    """Sample the Rust process's RSS and terminate it if the budget is exceeded."""
    started = time.monotonic()
    peak = 0
    stopped = False
    with log_path.open('w') as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=env)
        try:
            while process.poll() is None:
                try:
                    status = Path(f'/proc/{process.pid}/status').read_text()
                    fields = dict(line.split(':', 1) for line in status.splitlines() if ':' in line)
                    rss = int(fields.get('VmRSS', '0 kB').split()[0])
                    peak = max(peak, rss)
                    if rss > max_rss_mib * 1024:
                        stopped = True
                        process.terminate()
                        break
                except FileNotFoundError:
                    pass
                time.sleep(.02)
            try:
                code = process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill()
                code = process.wait()
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
    result = dict(exit_code=code, elapsed_seconds=round(time.monotonic()-started, 2),
                  peak_observed_rss_mib=round(peak/1024, 2), memory_limit_mib=max_rss_mib,
                  stopped_for_memory=stopped)
    atomic_json(log_path.with_suffix('.json'), result)
    if code:
        raise RuntimeError(f'{log_path.name}: {result}')
    return result


def grid_matches(source, other):
    return (source.shape == other.shape and source.crs == other.crs
            and source.transform.almost_equals(other.transform, precision=1e-10))


def compatible_tiff(folder):
    """Whitebox 2.4 cannot read floating-point TIFF predictor 3; preserve samples."""
    with rasterio.Env(GDAL_CACHEMAX=48*1024**2), rasterio.open(folder/'analysis.tif') as src:
        profile=src.profile.copy();profile.update(predictor=1)
        target=folder/'analysis-whitebox.tif'
        with rasterio.open(target,'w',**profile) as dst:
            for _,window in src.block_windows(1):
                dst.write(src.read(1,window=window),1,window=window)
        with rasterio.open(target) as dst:
            if not grid_matches(src,dst) or src.nodata!=dst.nodata:
                raise ValueError('Whitebox TIFF adaptation changed the grid or nodata')
            for _,window in src.block_windows(1):
                if not np.array_equal(src.read(1,window=window),dst.read(1,window=window),equal_nan=True):
                    raise ValueError('Whitebox TIFF adaptation changed elevation samples')


def installation(root, directory='tools/whitebox-prebuilt'):
    root=Path(root);record=read_json(root/directory/'build.json')
    binary=root/record['binary']
    if sha256(binary)!=record['binary_sha256']:raise ValueError('Whitebox binary changed')
    env=os.environ.copy()
    if record.get('library_directory'):
        libraries=root/record['library_directory']
        for name,digest in record['runtime_sha256'].items():
            if sha256(libraries/name)!=digest:raise ValueError('Whitebox runtime changed')
        env['LD_LIBRARY_PATH']=str(libraries)+(':'+env['LD_LIBRARY_PATH'] if env.get('LD_LIBRARY_PATH') else '')
    if record.get('settings') and read_json(binary.parent/'settings.json')!=record['settings']:
        raise ValueError('Whitebox settings changed')
    return record,binary,env


def run_flow(binary,folder,env,*,workflow='staged',logs=None,prefix='',max_rss_mib=1100):
    folder=Path(folder);logs=Path(logs or folder)
    if workflow=='full':
        commands=[('flow',[
            '-r=FlowAccumulationFullWorkflow','-v','--out_type=cells',
            f'--dem={folder/"analysis-whitebox.tif"}',f'--out_dem={folder/"conditioned.tif"}',
            f'--out_pntr={folder/"pointer.tif"}',f'--out_accum={folder/"accumulation.tif"}'])]
    elif workflow=='staged':
        commands=[
            ('breach',['-r=BreachDepressions','-v',f'--dem={folder/"analysis-whitebox.tif"}',
                       f'--output={folder/"conditioned.tif"}']),
            ('pointer',['-r=D8Pointer','-v',f'--dem={folder/"conditioned.tif"}',f'--output={folder/"pointer.tif"}']),
            ('accumulation',['-r=D8FlowAccumulation','-v','--pntr','--out_type=cells',
                             f'--input={folder/"pointer.tif"}',f'--output={folder/"accumulation.tif"}'])]
    else:raise ValueError('Unknown Whitebox workflow')
    for name,command in commands:
        yield name,measured([str(binary),*command],logs/(prefix+name+'.log'),max_rss_mib,env)


def route_context(root,folder,logs,prefix):
    """Export common compass directions and validity for the existing owner fit."""
    folder=Path(folder)
    _,binary,env=installation(root)
    compatible_tiff(folder)
    runs=dict(run_flow(binary,folder,env,logs=logs,prefix=prefix))
    with rasterio.open(folder/'analysis.tif') as src,rasterio.open(folder/'pointer.tif') as ds:
        if not grid_matches(src,ds):raise ValueError('Whitebox pointer grid changed')
        valid=src.read_masks(1)>0;values=ds.read(1)
        if not np.array_equal(ds.read_masks(1)>0,valid):raise ValueError('Whitebox pointer validity changed')
        if np.any(valid&~np.isin(values,[0,1,2,4,8,16,32,64,128])):
            raise ValueError('Unexpected Whitebox pointer code')
        pointer=np.where(valid,values,0).astype('uint8');del values
        uncertain=uncertain_upstream(pointer,valid)
        direction=TO_DIRECTION[pointer].astype('int16');direction[~valid]=-32768
        profile=src.profile.copy();profile.update(dtype='int16',nodata=-32768,predictor=2)
        with rasterio.open(folder/'drainage.tif','w',**profile) as dst:dst.write(direction,1)
        profile.update(dtype='uint8',nodata=None,predictor=1)
        with rasterio.open(folder/'evidence.tif') as evidence,rasterio.open(folder/'routing-valid.tif','w',**profile) as dst:
            if not grid_matches(src,evidence):raise ValueError('Evidence grid changed')
            dst.write(((evidence.read(1)>0)&valid&~uncertain).astype('uint8'),1)
        # Every routing output must preserve the original analysis NULL footprint.
        for name in ('conditioned.tif','accumulation.tif'):
            with rasterio.open(folder/name) as other:
                if not grid_matches(src,other):raise ValueError(name+': grid changed')
                for _,window in src.block_windows(1):
                    if not np.array_equal(src.read_masks(1,window=window)>0,other.read_masks(1,window=window)>0):
                        raise ValueError(name+': nodata footprint changed')
    return dict(runs=runs,uncertain_analysis_pixels=int(uncertain.sum()))

