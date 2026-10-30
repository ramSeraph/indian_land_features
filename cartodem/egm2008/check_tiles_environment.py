"""Verify the small tiling runtime, including GDAL, lossless WebP and opacity."""
import argparse
import importlib.metadata
import json
import platform
from pathlib import Path

import numpy as np
import rasterio
from PIL import features
from pmtiles.tile import tileid_to_zxy,zxy_to_tileid
from production_tiles import validate_rgba,webp_bytes
from terrain_rgb import encode_terrain_rgb


def check():
    if not features.check('webp'):
        raise RuntimeError('Pillow lacks WebP support')
    with rasterio.Env() as env:
        for driver in ('GTiff','COG','VRT'):
            if driver not in env.drivers():raise RuntimeError(f'GDAL lacks {driver}')
    heights=np.array([[0,16,np.nan,-255]],dtype='float32')
    valid=np.isfinite(heights)
    rgba=encode_terrain_rgb(heights,round_bits=4)
    validate_rgba(rgba,12,valid)
    np.testing.assert_array_equal(rgba[0,0],rgba[0,2])
    np.testing.assert_array_equal(rgba[0,2],[1,134,160,255])
    webp_bytes(rgba)  # Also decodes and compares every channel.
    assert tileid_to_zxy(zxy_to_tileid(12,2902,1908))==(12,2902,1908)
    versions={p:importlib.metadata.version(p) for p in
              ('numpy','rasterio','pillow','pmtiles','affine','attrs','certifi','click','pyparsing')}
    return dict(python=platform.python_version(),system=platform.system(),
                architecture=platform.machine(),gdal=rasterio.__gdal_version__,
                proj=rasterio.__proj_version__,packages=versions,
                nodata_rgba=[1,134,160,255],lossless_webp_roundtrip=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=check()
    data=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(data)
    print(data,end='')
