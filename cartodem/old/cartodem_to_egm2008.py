#!/usr/bin/env python3
"""Convert original CartoDEM Int16 sheets from assumed ellipsoidal heights to EGM2008.

Usage:
  python cartodem_to_egm2008.py INPUT.tif GEOID.tif OUTPUT.tif
Uses only CartoDEM and the pinned EGM2008 grid. No comparison or boundary mask.
"""
import os
os.environ.setdefault("GDAL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import hashlib
import json
import tempfile
from pathlib import Path
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window

GEOID_SHA256 = "4191d471eefebf24091b56dbc604353cb3b8cf8cc70e448bb9ae56a272bef17a"

def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def windows(width, height):
    for y in range(0, height, 512):
        for x in range(0, width, 512):
            yield Window(x, y, min(512, width-x), min(512, height-y))

def convert(source_path, geoid_path, output_path, source_nodata=-255):
    source_path, geoid_path, output_path = map(Path, (source_path, geoid_path, output_path))
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite {output_path}")
    if output_path.resolve() in (source_path.resolve(), geoid_path.resolve()):
        raise ValueError("Output must differ from both inputs")
    geoid_hash = sha256(geoid_path)
    if geoid_hash != GEOID_SHA256:
        raise ValueError("Geoid differs from the pinned us_nga_egm08_25.tif; verify its model, units and sign")
    source_hash = sha256(source_path)
    temporary = None
    valid_count = nodata_count = geoid_missing = 0
    nmin, nmax = np.inf, -np.inf
    hmin, hmax = np.inf, -np.inf
    try:
        with rasterio.Env(GDAL_CACHEMAX=24*1024**2, GDAL_NUM_THREADS="1"), \
             rasterio.open(source_path) as src, rasterio.open(geoid_path) as geo:
            if src.count != 1 or src.dtypes != ("int16",) or src.crs != rasterio.crs.CRS.from_epsg(4326):
                raise ValueError("Expected an original, single-band Int16 CartoDEM sheet in EPSG:4326")
            if src.scales != (1.0,) or src.offsets != (0.0,):
                raise ValueError("Unexpected scale/offset; inspect the input before conversion")
            if src.transform.b != 0 or src.transform.d != 0:
                raise ValueError("Expected a north-up original sheet")
            if "HEIGHT_FORMULA" in src.tags() or "VERTICAL_DATUM" in src.tags():
                raise ValueError("Input already has height-datum processing tags; inspect before converting")
            profile = dict(driver="GTiff", width=src.width, height=src.height, count=1,
                           dtype="float32", crs=src.crs, transform=src.transform,
                           nodata=np.nan, tiled=True, blockxsize=256, blockysize=256,
                           compress="DEFLATE", predictor=3, BIGTIFF="IF_SAFER")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(prefix=output_path.stem+".", suffix=".part.tif",
                                        dir=output_path.parent)
            os.close(fd)
            temporary = Path(name)
            # Sample scalar N only. The grid's EPSG:4979 metadata must not
            # trigger an automatic vertical transform during this resampling.
            with WarpedVRT(geo, src_crs="EPSG:4326", crs=src.crs,
                           transform=src.transform, width=src.width, height=src.height,
                           dtype="float32", resampling=Resampling.bilinear,
                           nodata=np.nan, warp_mem_limit=16, tolerance=1e-6) as geoid, \
                 rasterio.open(temporary, "w", **profile) as dst:
                dst.set_band_description(1, "EGM2008 height under WGS84 ellipsoidal source assumption")
                dst.set_band_unit(1, "m")
                if "AREA_OR_POINT" in src.tags():
                    dst.update_tags(AREA_OR_POINT=src.tags()["AREA_OR_POINT"])
                dst.update_tags(
                    HEIGHT_FORMULA="H = h - N",
                    SOURCE_VERTICAL_DATUM="Assumed WGS84 ellipsoidal heights in metres",
                    VERTICAL_DATUM="EGM2008 under stated source assumption",
                    TARGET_VERTICAL_CRS="EPSG:3855",
                    GEOID_MODEL="EGM2008",
                    GEOID_GRID="us_nga_egm08_25.tif",
                    GEOID_SHA256=geoid_hash,
                    GEOID_RESAMPLING="bilinear at original CartoDEM pixel centres",
                    SOURCE_FILE=source_path.name,
                    SOURCE_SHA256=source_hash,
                    SOURCE_NODATA_POLICY=f"existing mask or nonfinite or original value == {source_nodata}",
                    RASTERIO_VERSION=rasterio.__version__,
                    GDAL_VERSION=rasterio.__gdal_version__,
                    PROJ_VERSION=rasterio.__proj_version__,
                )
                for window in windows(src.width, src.height):
                    h = src.read(1, window=window).astype(np.float32)
                    n = geoid.read(1, window=window)
                    source_valid = (src.read_masks(1, window=window) != 0) & np.isfinite(h) & (h != source_nodata)
                    valid = source_valid & np.isfinite(n)
                    result = np.full(h.shape, np.nan, dtype=np.float32)
                    np.subtract(h, n, out=result, where=valid)
                    dst.write(result, 1, window=window)
                    valid_count += int(valid.sum())
                    nodata_count += int((~source_valid).sum())
                    geoid_missing += int((source_valid & ~np.isfinite(n)).sum())
                    if valid.any():
                        nmin = min(nmin, float(n[valid].min()))
                        nmax = max(nmax, float(n[valid].max()))
                        hmin = min(hmin, float(result[valid].min()))
                        hmax = max(hmax, float(result[valid].max()))
            if geoid_missing:
                raise ValueError(f"Geoid unavailable at {geoid_missing} otherwise-valid input pixels")
            with rasterio.open(temporary) as check:
                assert check.transform == src.transform and check.crs == src.crs and check.shape == src.shape
                assert check.dtypes == ("float32",) and np.isnan(check.nodata)
                written_valid = sum(int(np.isfinite(check.read(1, window=w)).sum())
                                    for w in windows(check.width, check.height))
                assert written_valid == valid_count
            if output_path.exists():
                raise FileExistsError(f"Output appeared during conversion: {output_path}")
            temporary.replace(output_path)
            temporary = None
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return dict(source=str(source_path), output=str(output_path),
                source_sha256=source_hash, geoid_sha256=geoid_hash,
                valid_pixels=valid_count, source_nodata_pixels=nodata_count,
                geoid_missing_pixels=geoid_missing,
                geoid_undulation_range_m=[nmin,nmax] if valid_count else None,
                corrected_height_range_m=[hmin,hmax] if valid_count else None,
                formula="H = h - N",
                vertical_assumption="Original CartoDEM is WGS84 ellipsoidal; target EGM2008",
                versions=dict(numpy=np.__version__, rasterio=rasterio.__version__,
                              gdal=rasterio.__gdal_version__, proj=rasterio.__proj_version__))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("geoid", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--source-nodata", type=float, default=-255,
                        help="Inferred fill in this original CartoDEM release (default: -255)")
    args = parser.parse_args()
    print(json.dumps(convert(args.source, args.geoid, args.output, args.source_nodata),
                     indent=2, allow_nan=False))

if __name__ == "__main__":
    main()
