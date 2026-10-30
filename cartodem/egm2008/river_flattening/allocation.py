"""Assign river cells to centreline seeds using Whitebox raster cost allocation.

Only disposable rasters use a metric plane. Published grids and elevations never
pass through Whitebox here. Unit cost, eight neighbours and the existing dx/dy
preserve the old graph's path metric, including barriers around land and islands.
"""
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

from build import atomic_json, sha256
from .whitebox import grid_matches, installation, measured


def backend_recipe(root, backend):
    if backend not in ('whitebox', 'scipy'):
        raise ValueError('Unknown river allocation backend')
    result = dict(backend=backend, code_sha256=sha256(Path(__file__)))
    if backend == 'whitebox':
        result.update(tools=['CostDistance', 'CostAllocation'],
                      whitebox=installation(root)[0],
                      adapter_sha256=sha256(Path(__file__).with_name('whitebox.py')),
                      neighbours=8, cost=1, outside='nodata barrier',
                      metric='existing constant window dx/dy in metres')
    return result


def write_inputs(folder, river, seeds, dx, dy):
    """Write unit cost and unique positive source IDs, in bounded TIFF blocks."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    if (river.ndim != 2 or seeds.shape != river.shape or river.dtype != bool
            or seeds.dtype != bool or not np.isfinite([dx, dy]).all() or min(dx, dy) <= 0):
        raise ValueError('Invalid allocation inputs')
    seed_pixels = np.flatnonzero(seeds)
    if not len(seed_pixels) or np.any(seeds & ~river):
        raise ValueError('Invalid centreline seeds')
    if len(seed_pixels) > np.iinfo('int32').max:
        raise ValueError('Too many centreline seeds for Int32 labels')
    height, width = river.shape
    profile = dict(driver='GTiff', width=width, height=height, count=1,
                   crs='EPSG:3857', transform=from_origin(0, height*dy, dx, dy),
                   tiled=True, blockxsize=256, blockysize=256,
                   compress='DEFLATE', predictor=1)
    with rasterio.open(folder/'cost.tif', 'w', **profile, dtype='uint8', nodata=255) as cost, \
         rasterio.open(folder/'sources.tif', 'w', **profile, dtype='int32', nodata=-1) as sources:
        for _, window in cost.block_windows(1):
            y, x, h, w = map(int, (window.row_off, window.col_off, window.height, window.width))
            water = river[y:y+h, x:x+w]
            cost.write(np.where(water, 1, 255).astype('uint8'), 1, window=window)
            labels = np.where(water, 0, -1).astype('int32')
            r, c = np.nonzero(seeds[y:y+h, x:x+w])
            labels[r, c] = np.searchsorted(seed_pixels, (r+y)*width+c+x)+1
            sources.write(labels, 1, window=window)
        sources.update_tags(COORDINATES='Disposable metric plane; not geographic output',
                            SEED_LABELS='1-based index into row-major centreline pixels')
    return dict(shape=list(river.shape), dx_m=float(dx), dy_m=float(dy),
                river_pixels=int(river.sum()), seeds=len(seed_pixels),
                files_sha256={p:sha256(folder/p) for p in ('cost.tif', 'sources.tif')})


def run_tools(root, folder):
    tool, binary, env = installation(root)
    results = {}
    for name, args in (
        ('cost-distance', ['-r=CostDistance', f'--source={folder/"sources.tif"}',
                           f'--cost={folder/"cost.tif"}', f'--out_accum={folder/"distance.tif"}',
                           f'--out_backlink={folder/"backlink.tif"}']),
        ('cost-allocation', ['-r=CostAllocation', f'--source={folder/"sources.tif"}',
                             f'--backlink={folder/"backlink.tif"}', f'--output={folder/"owners.tif"}'])):
        results[name] = measured([str(binary), '-v', *args], folder/(name+'.log'), 1100, env)
    return tool, results


def read_owners(folder, river, seeds):
    """Validate coverage/seed identity, returning one local seed pixel per river cell."""
    flat = np.flatnonzero(river)
    seed_pixels = np.flatnonzero(seeds)
    owners = np.empty(len(flat), dtype='int64')
    width = river.shape[1]
    with rasterio.open(folder/'cost.tif') as cost, \
         rasterio.open(folder/'owners.tif') as output:
        if not grid_matches(cost, output) or output.shape != river.shape:
            raise ValueError('Allocation grid changed')
        for _, window in cost.block_windows(1):
            y, x, h, w = map(int, (window.row_off, window.col_off, window.height, window.width))
            water = river[y:y+h, x:x+w]
            labels = output.read(1, window=window)
            values = labels[water]
            if (np.any(~np.isfinite(values)) or np.any(values != np.floor(values))
                    or np.any((values < 1) | (values > len(seed_pixels)))):
                raise ValueError('A river pixel has no valid centreline owner')
            if not np.array_equal(output.read_masks(1, window=window) > 0, water):
                raise ValueError('Allocation crossed the river mask')
            r, c = np.nonzero(seeds[y:y+h, x:x+w])
            if not np.array_equal(labels[r, c], np.searchsorted(seed_pixels, (r+y)*width+c+x)+1):
                raise ValueError('A centreline seed was assigned to another seed')
            r, c = np.nonzero(water)
            owners[np.searchsorted(flat, (r+y)*width+c+x)] = seed_pixels[values.astype('int64')-1]
    return owners


def make_allocator(root, work, sheet, recipe):
    """Small reference tests can explicitly request the historical SciPy backend."""
    backend = recipe['allocation']['backend']
    if backend == 'scipy':
        return None
    if backend != 'whitebox':
        raise ValueError('Unknown river allocation backend')

    def allocate(river, seeds, dx, dy):
        folder = work/'allocation'/Path(sheet).stem
        # Retry scratch outputs from the same recipe by recomputing both tools.
        # No allocation result is reused without the geometry checkpoint check.
        info = write_inputs(folder, river, seeds, dx, dy)
        atomic_json(folder/'inputs.json', info)
        tool, runs = run_tools(root, folder)
        if tool != recipe['allocation']['whitebox']:
            raise ValueError('Whitebox installation changed after geometry setup')
        owners = read_owners(folder, river, seeds)
        atomic_json(folder/'report.json', dict(
            state='complete', inputs=info, runs=runs, assigned_river_pixels=len(owners),
            allocation=recipe['allocation'],
            outputs_sha256={p:sha256(folder/p) for p in ('distance.tif', 'backlink.tif', 'owners.tif')},
            mask_preserved=True, seed_identity_preserved=True))
        # Only this adapter's disposable rasters are removed; provenance/logs stay.
        for name in ('cost.tif', 'sources.tif', 'distance.tif', 'backlink.tif', 'owners.tif'):
            (folder/name).unlink()
        return owners

    return allocate
