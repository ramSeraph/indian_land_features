#!/usr/bin/env python3
"""Build a Terrain-RGBA XYZ pyramid from a corrected EPSG:4326 elevation raster.

Use a mosaic of corrected sheets for production. A single sheet is useful for a
pilot, but its edge tiles contain nodata where neighboring sheets are absent.
"""
import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.warp import transform, transform_bounds
from rasterio.windows import Window

from terrain_rgb import encode_terrain_rgb, write_terrain_rgb

HALF_WORLD = math.pi * 6378137


def tile_pixel(lon, lat, zoom, tile_size):
    tx = (lon + 180) / 360 * 2**zoom
    ty = (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * 2**zoom
    return math.floor(tx), math.floor(ty), int(tx % 1 * tile_size), int(ty % 1 * tile_size)


def build_tiles(source_path, output_dir, min_zoom=5, max_zoom=12, tile_size=512,
                sample=None, precision="full"):
    source_path, output_dir = Path(source_path), Path(output_dir)
    if not 0 <= min_zoom <= max_zoom <= 24:
        raise ValueError("Expected 0 <= min_zoom <= max_zoom <= 24")
    if tile_size not in (256, 512):
        raise ValueError("Tile size must be 256 or 512")
    if precision not in ("full", "legacy"):
        raise ValueError("precision must be full or legacy")
    if precision == "legacy" and not 5 <= min_zoom <= max_zoom <= 12:
        raise ValueError("The legacy precision schedule is defined for zooms 5 through 12")
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    started = time.perf_counter()
    report = {"source": str(source_path), "output": str(output_dir),
              "tile_size": tile_size, "minzoom": min_zoom, "maxzoom": max_zoom,
              "encoding": "mapbox", "format": "webp", "alpha": 255,
              "nodata_preserved": False, "nodata_rgb": [1, 134, 160],
              "precision": precision, "zooms": []}
    with rasterio.Env(GDAL_CACHEMAX=24*1024**2, GDAL_NUM_THREADS="1"), \
         rasterio.open(source_path) as src:
        if (src.count != 1 or src.dtypes != ("float32",)
                or src.crs != rasterio.crs.CRS.from_epsg(4326)
                or src.nodata is None or not np.isnan(src.nodata)):
            raise ValueError("Expected corrected single-band Float32 EPSG:4326 with NaN nodata")
        bounds = transform_bounds(src.crs, "EPSG:3857", *src.bounds, densify_pts=21)
        if not (-HALF_WORLD <= bounds[1] < bounds[3] <= HALF_WORLD):
            raise ValueError("Raster exceeds Web Mercator latitude limits")
        cx, cy = src.xy(src.height // 2, src.width // 2)
        px, py = transform(src.crs, "EPSG:3857",
                           [cx, cx + src.transform.a, cx],
                           [cy, cy, cy + src.transform.e])
        native_resolution = max(abs(px[1] - px[0]), abs(py[2] - py[0]))
        report["bounds"] = list(src.bounds)
        report["native_mercator_pixel_size_m_at_center"] = native_resolution
        output_dir.mkdir(parents=True)
        for z in range(min_zoom, max_zoom + 1):
            zoom_started = time.perf_counter()
            round_bits = 16 - z if precision == "legacy" else 0
            span = 2 * HALF_WORLD / 2**z
            resolution = span / tile_size
            xmin = max(0, math.floor((bounds[0] + HALF_WORLD) / span))
            xmax = min(2**z - 1, math.ceil((bounds[2] + HALF_WORLD) / span) - 1)
            ymin = max(0, math.floor((HALF_WORLD - bounds[3]) / span))
            ymax = min(2**z - 1, math.ceil((HALF_WORLD - bounds[1]) / span) - 1)
            method = Resampling.average if resolution >= native_resolution else Resampling.bilinear
            stats = {"zoom": z, "resampling": method.name, "tiles": 0,
                     "bytes": 0, "valid_pixels": 0, "nodata_pixels": 0,
                     "fully_nodata_tiles": 0, "round_bits": round_bits,
                     "height_step_m": 0.1 * 2**round_bits,
                     "max_encoding_error_m": 0.0}
            wanted = tile_pixel(*sample, z, tile_size) if sample else None
            # A common grid for every tile at this zoom prevents per-tile alignment drift.
            # Explicit 2-D CRSs keep this a horizontal reprojection only.
            with WarpedVRT(src, src_crs="EPSG:4326", crs="EPSG:3857",
                           transform=from_origin(-HALF_WORLD + xmin*span,
                                                 HALF_WORLD - ymin*span,
                                                 resolution, resolution),
                           width=(xmax-xmin+1)*tile_size, height=(ymax-ymin+1)*tile_size,
                           dtype="float32", src_nodata=np.nan, nodata=np.nan,
                           resampling=method, warp_mem_limit=16, tolerance=1e-6) as vrt:
                for x in range(xmin, xmax + 1):
                    for y in range(ymin, ymax + 1):
                        window = Window((x-xmin)*tile_size, (y-ymin)*tile_size,
                                        tile_size, tile_size)
                        heights = vrt.read(1, window=window, masked=True)
                        valid = ~np.ma.getmaskarray(heights) & np.isfinite(heights.data)
                        count = int(valid.sum())
                        path = output_dir / str(z) / str(x) / f"{y}.webp"
                        write_terrain_rgb(heights, path, round_bits=round_bits)
                        rgba = encode_terrain_rgb(heights, round_bits=round_bits)
                        rgb = rgba[..., :3].astype(np.int32)
                        decoded_heights = -10000 + .1 * (rgb[..., 0]*65536
                                                        + rgb[..., 1]*256 + rgb[..., 2])
                        if count:
                            error = float(np.max(np.abs(decoded_heights[valid] - heights.data[valid])))
                            if error > stats["height_step_m"] / 2 + 1e-6:
                                raise ValueError("Encoding exceeded the selected precision's error bound")
                            stats["max_encoding_error_m"] = max(stats["max_encoding_error_m"], error)
                        stats["tiles"] += 1
                        stats["bytes"] += path.stat().st_size
                        stats["valid_pixels"] += count
                        stats["nodata_pixels"] += tile_size**2 - count
                        stats["fully_nodata_tiles"] += int(count == 0)
                        if wanted and (x, y) == wanted[:2]:
                            col, row = wanted[2:]
                            pixel_rgba = rgba[row, col].tolist()
                            decoded = float(decoded_heights[row, col])
                            numeric = float(heights.data[row, col]) if valid[row, col] else None
                            stats["sample"] = {"tile": [z, x, y], "pixel": [col, row],
                                               "rgba": pixel_rgba, "resampled_height_m": numeric,
                                               "decoded_height_m": decoded, "source_valid": bool(valid[row,col])}
            stats["elapsed_seconds"] = time.perf_counter() - zoom_started
            report["zooms"].append(stats)
            print(json.dumps(stats, allow_nan=False), flush=True)
    report["tiles"] = sum(z["tiles"] for z in report["zooms"])
    report["image_bytes"] = sum(z["bytes"] for z in report["zooms"])
    report["elapsed_seconds"] = time.perf_counter() - started
    metadata = {"tilejson": "3.0.0", "name": source_path.stem, "scheme": "xyz",
                "tiles": ["{z}/{x}/{y}.webp"], "format": "webp", "encoding": "mapbox",
                "tileSize": tile_size, "minzoom": min_zoom, "maxzoom": max_zoom,
                "bounds": report["bounds"], "alpha": 255, "nodata_height_m": 0, "nodata_preserved": False}
    (output_dir / "tiles.json").write_text(json.dumps(metadata, indent=2))
    (output_dir / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--min-zoom", type=int, default=5)
    parser.add_argument("--max-zoom", type=int, default=12)
    parser.add_argument("--tile-size", type=int, choices=(256, 512), default=512)
    parser.add_argument("--sample", type=float, nargs=2, metavar=("LON", "LAT"))
    parser.add_argument("--precision", choices=("full", "legacy"), default="full",
                        help="full: 0.1 m; legacy: erase 11 bits at z5 through 4 bits at z12")
    args = parser.parse_args()
    report = build_tiles(args.source, args.output, args.min_zoom, args.max_zoom,
                         args.tile_size, args.sample, args.precision)
    print(json.dumps({k: v for k, v in report.items() if k != "zooms"}, indent=2))


if __name__ == "__main__":
    main()
