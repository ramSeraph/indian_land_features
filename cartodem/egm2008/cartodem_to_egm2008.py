#!/usr/bin/env python3
"""Convert original CartoDEM Int16 sheets from assumed ellipsoidal heights to EGM2008.

Usage:
  python cartodem_to_egm2008.py INPUT.tif GEOID.tif OUTPUT.tif
Both original fill values are masked before correction. Optional GLO-30 WBM
ocean masking and finite output nodata support the production coastal build.
"""
import os
os.environ.setdefault("GDAL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import hashlib
import json
import tempfile
from contextlib import ExitStack
from pathlib import Path
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window

GEOID_SHA256 = "4191d471eefebf24091b56dbc604353cb3b8cf8cc70e448bb9ae56a272bef17a"
SOURCE_NODATA_VALUES = (-255, -1024)

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

def convert(source_path, geoid_path, output_path, source_nodata=SOURCE_NODATA_VALUES,
            *, water_mask_path=None, output_nodata=np.nan):
    source_path, geoid_path, output_path = map(Path, (source_path, geoid_path, output_path))
    source_nodata = np.atleast_1d(source_nodata).astype(np.float64).tolist()
    if not source_nodata or not np.isfinite(source_nodata).all():
        raise ValueError("Source nodata values must be a nonempty list of finite numbers")
    output_nodata = float(np.float32(output_nodata))
    if np.isinf(output_nodata):
        raise ValueError("Output nodata must be NaN or a finite Float32 value")
    water_mask_path = Path(water_mask_path) if water_mask_path is not None else None
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite {output_path}")
    inputs = [source_path, geoid_path] + ([water_mask_path] if water_mask_path else [])
    if output_path.resolve() in [p.resolve() for p in inputs]:
        raise ValueError("Output must differ from both inputs")
    geoid_hash = sha256(geoid_path)
    if geoid_hash != GEOID_SHA256:
        raise ValueError("Geoid differs from the pinned us_nga_egm08_25.tif; verify its model, units and sign")
    source_hash = sha256(source_path)
    mask_hash = sha256(water_mask_path) if water_mask_path else None
    temporary = None
    valid_count = nodata_count = geoid_missing = 0
    source_valid_count = ocean_count = removed_by_ocean = missing_on_land = 0
    source_fill_counts = {str(int(v)): 0 for v in source_nodata}
    water_counts = {str(v): 0 for v in range(4)}
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
                           nodata=output_nodata, tiled=True, blockxsize=256, blockysize=256,
                           compress="DEFLATE", predictor=3, BIGTIFF="IF_SAFER")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(prefix=output_path.stem+".", suffix=".part.tif",
                                        dir=output_path.parent)
            os.close(fd)
            temporary = Path(name)
            # Sample scalar N only. The grid's EPSG:4979 metadata must not
            # trigger an automatic vertical transform during this resampling.
            with ExitStack() as stack:
                geoid = stack.enter_context(WarpedVRT(geo, src_crs="EPSG:4326", crs=src.crs,
                           transform=src.transform, width=src.width, height=src.height,
                           dtype="float32", resampling=Resampling.bilinear,
                           nodata=np.nan, warp_mem_limit=16, tolerance=1e-6))
                water = None
                if water_mask_path:
                    mask_src = stack.enter_context(rasterio.open(water_mask_path))
                    if mask_src.count != 1 or mask_src.crs is None:
                        raise ValueError("Expected a georeferenced single-band water mask")
                    # WBM class zero means land, even when its TIFF tags it nodata.
                    water = stack.enter_context(WarpedVRT(
                        mask_src, crs=src.crs, transform=src.transform,
                        width=src.width, height=src.height, dtype="uint8",
                        src_nodata=255, nodata=255, resampling=Resampling.nearest,
                        warp_mem_limit=16, tolerance=1e-6))
                dst = stack.enter_context(rasterio.open(temporary, "w", **profile))
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
                    SOURCE_NODATA_POLICY=f"existing mask or nonfinite or original value in {source_nodata}",
                    SOURCE_NODATA_VALUES=json.dumps(source_nodata),
                    RASTERIO_VERSION=rasterio.__version__,
                    GDAL_VERSION=rasterio.__gdal_version__,
                    PROJ_VERSION=rasterio.__proj_version__,
                    OUTPUT_NODATA=str(output_nodata),
                    COASTAL_MASK="GLO30 WBM ocean class 1" if water else "none",
                    WATER_MASK_FILE=water_mask_path.name if water_mask_path else "",
                    WATER_MASK_SHA256=mask_hash or "",
                )
                for window in windows(src.width, src.height):
                    h = src.read(1, window=window).astype(np.float32)
                    n = geoid.read(1, window=window)
                    source_valid = ((src.read_masks(1, window=window) != 0)
                                    & np.isfinite(h) & ~np.isin(h, source_nodata))
                    valid = source_valid & np.isfinite(n)
                    source_valid_count += int(source_valid.sum())
                    for v in source_nodata:
                        source_fill_counts[str(int(v))] += int((h == v).sum())
                    if water:
                        classes = water.read(1, window=window, masked=False)
                        if np.any(classes > 3):
                            raise ValueError("Water mask has missing coverage or unknown classes")
                        ocean = classes == 1
                        ocean_count += int(ocean.sum())
                        removed_by_ocean += int((source_valid & ocean).sum())
                        missing_on_land += int((~source_valid & (classes == 0)).sum())
                        for v in range(4):
                            water_counts[str(v)] += int((classes == v).sum())
                        valid &= ~ocean
                    result = np.full(h.shape, output_nodata, dtype=np.float32)
                    np.subtract(h, n, out=result, where=valid)
                    if np.isfinite(output_nodata) and np.any(valid & (result == output_nodata)):
                        raise ValueError(f"Valid corrected height collides with output nodata {output_nodata}")
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
                assert check.dtypes == ("float32",)
                assert ((np.isnan(output_nodata) and np.isnan(check.nodata))
                        or check.nodata == output_nodata)
                written_valid = sum(int((check.read_masks(1, window=w) != 0).sum())
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
                source_valid_pixels=source_valid_count,
                output_nodata_pixels=nodata_count + removed_by_ocean,
                source_fill_counts=source_fill_counts,
                water_mask=str(water_mask_path) if water_mask_path else None,
                water_mask_sha256=mask_hash,
                water_class_counts=water_counts if water_mask_path else None,
                ocean_pixels=ocean_count, newly_masked_ocean_pixels=removed_by_ocean,
                source_nodata_on_glo30_land_pixels=missing_on_land,
                output_nodata=output_nodata if np.isfinite(output_nodata) else "NaN",
                source_nodata_values=source_nodata,
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
    parser.add_argument("--source-nodata", type=float, nargs="+", default=SOURCE_NODATA_VALUES,
                        help="Original fill values to mask before correction (default: -255 -1024)")
    parser.add_argument("--water-mask", type=Path, help="GLO-30 WBM; class 1 becomes nodata")
    parser.add_argument("--output-nodata", type=float, default=np.nan,
                        help="Output sentinel, default NaN; production uses -255")
    args = parser.parse_args()
    print(json.dumps(convert(args.source, args.geoid, args.output, args.source_nodata,
                             water_mask_path=args.water_mask, output_nodata=args.output_nodata),
                     indent=2, allow_nan=False))

if __name__ == "__main__":
    main()
