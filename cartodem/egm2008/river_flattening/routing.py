"""Prepare temporary routing heights while keeping original nodata out of evidence."""
import numpy as np
from rasterio.fill import fillnodata
from scipy.ndimage import binary_dilation

# Normalized D8 codes: NE, N, NW, W, SW, S, SE, E.
OFFSETS = np.array([[0,0],[-1,1],[-1,0],[-1,-1],[0,-1],
                    [1,-1],[1,0],[1,1],[0,1]],dtype=np.int8)


def analysis_window(height, classes, max_search=64):
    valid = np.isfinite(height)&(height!=-255)&(classes!=1)&(classes<4)
    h = np.where(valid,height,0).astype(np.float32)
    if valid.any():
        h = fillnodata(h, mask=valid.astype('uint8'), max_search_distance=max_search)
        # GDAL leaves distant holes at zero: identify fillable cells separately.
        from scipy.ndimage import distance_transform_edt
        filled = (~valid)&(classes!=1)&(classes<4)&(distance_transform_edt(~valid)<=max_search)
    else:
        filled = np.zeros(valid.shape,bool)
    h[~(valid|filled)] = -255
    # Exclude water-flow samples near either provisional interpolation or holes.
    uncertain = binary_dilation((~valid)&(classes!=1), iterations=8)
    evidence = valid&~uncertain
    return h,evidence,filled
