"""Export final native-grid Float32 elevations as lossless COGs, verifying pixels."""
from pathlib import Path

import numpy as np
import rasterio
from rasterio.shutil import copy as raster_copy

from cartodem_to_egm2008 import sha256,windows


def export_cog(source,output):
    source,output=Path(source),Path(output)
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    temp=output.with_name(output.name+'.part')
    try:
        with rasterio.Env(GDAL_CACHEMAX=64*1024**2,GDAL_NUM_THREADS='1'):
            with rasterio.open(source) as src:
                if src.count!=1 or src.dtypes!=('float32',) or src.nodata!=-255:
                    raise ValueError('Final elevation sheets must be Float32 with NoData=-255')
            raster_copy(source,temp,driver='COG',COMPRESS='DEFLATE',PREDICTOR='FLOATING_POINT',
                        BLOCKSIZE=512,BIGTIFF='IF_SAFER',NUM_THREADS='1',
                        OVERVIEWS='IGNORE_EXISTING',OVERVIEW_RESAMPLING='AVERAGE')
            valid_count=0
            with rasterio.open(source) as src,rasterio.open(temp) as dst:
                if (src.shape!=dst.shape or src.transform!=dst.transform or src.crs!=dst.crs
                        or dst.nodata!=-255 or dst.dtypes!=('float32',)
                        or dst.tags(ns='IMAGE_STRUCTURE').get('LAYOUT')!='COG'):
                    raise ValueError('COG geometry, nodata or layout differs from the prepared sheet')
                if max(src.shape)>512 and not dst.overviews(1):
                    raise ValueError('Expected numeric COG overviews')
                for window in windows(src.width,src.height):
                    original=src.read(1,window=window)
                    result=dst.read(1,window=window)
                    if not np.isfinite(result).all() or not np.array_equal(original,result):
                        raise ValueError('COG changed native elevations or contains nonfinite samples')
                    original_mask=src.read_masks(1,window=window)
                    if not np.array_equal(original_mask,dst.read_masks(1,window=window)):
                        raise ValueError('COG changed native validity')
                    if np.any(result[original_mask==0]!=-255):
                        raise ValueError('COG nodata pixels were not uniformly stored as -255')
                    valid_count+=int(np.sum(original_mask!=0))
                overviews=dst.overviews(1)
        temp.replace(output)
    finally:
        temp.unlink(missing_ok=True)
    return dict(file=output.name,bytes=output.stat().st_size,sha256=sha256(output),
                valid_pixels=valid_count,nodata=-255,overviews=overviews,
                native_pixels_verified=True,layout='COG',overview_resampling='AVERAGE')
