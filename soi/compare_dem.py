# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "rasterio>=1.4", "scipy>=1.14", "pyproj>=3.6", "matplotlib>=3.9"]
# ///
"""Compare an SOI DEM with a window of public Copernicus GLO-30.

Numeric heights are preserved: this does not infer or convert SOI's vertical datum.
Remote COG reads fetch the sheet window, not the full global tileset.
"""

import argparse
import csv
import hashlib
import json
import math
from contextlib import ExitStack
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.merge import merge
from rasterio.warp import Resampling, reproject, transform_bounds
from rasterio.windows import Window, bounds as window_bounds, from_bounds
from scipy.ndimage import binary_erosion, map_coordinates


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def stats(values):
    values = np.asarray(values, dtype="float64")
    values = values[np.isfinite(values)]
    if not values.size:
        return {"count": 0}
    median = float(np.median(values))
    return {
        "count": int(values.size),
        "mean_m": float(values.mean()),
        "median_m": median,
        "mae_m": float(np.abs(values).mean()),
        "rmse_m": float(np.sqrt(np.mean(values**2))),
        "std_m": float(values.std()),
        "nmad_m": float(1.4826 * np.median(np.abs(values - median))),
        "p05_m": float(np.percentile(values, 5)),
        "p95_m": float(np.percentile(values, 95)),
        "min_m": float(values.min()),
        "max_m": float(values.max()),
    }


def sample_bilinear(data, transform, x, y):
    # Raster transforms locate pixel corners; array indices locate pixel centres.
    col, row = ~transform * (np.asarray(x), np.asarray(y))
    return map_coordinates(
        data,
        [row - 0.5, col - 0.5],
        order=1,
        mode="constant",
        cval=np.nan,
    )


def align_bounds(bounds, transform):
    """Enclose requested bounds without moving the native raster's pixel grid."""
    window = from_bounds(*bounds, transform=transform)
    left, top = math.floor(window.col_off), math.floor(window.row_off)
    right = math.ceil(window.col_off + window.width)
    bottom = math.ceil(window.row_off + window.height)
    return window_bounds(Window(left, top, right - left, bottom - top), transform)


def write_grid(path, data, profile, description, **tags):
    profile = dict(profile)
    profile.update(
        driver="GTiff",
        count=1,
        dtype="float32",
        nodata=-9999,
        compress="deflate",
        tiled=True,
        blockxsize=256,
        blockysize=256,
    )
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(np.where(np.isfinite(data), data, -9999).astype("float32"), 1)
        dst.set_band_description(1, description)
        dst.update_tags(**tags)


def reference_window(dem, output):
    # Allow enough margin for reprojection and bilinear sampling at sheet edges.
    b = dem.bounds
    bounds = transform_bounds(
        dem.crs,
        "EPSG:4326",
        b.left - 150,
        b.bottom - 150,
        b.right + 150,
        b.top + 150,
        densify_pts=41,
    )
    urls = []
    for lat in range(math.floor(bounds[1]), math.floor(bounds[3]) + 1):
        for lon in range(math.floor(bounds[0]), math.floor(bounds[2]) + 1):
            ns, ew = ("N" if lat >= 0 else "S"), ("E" if lon >= 0 else "W")
            name = (
                f"Copernicus_DSM_COG_10_{ns}{abs(lat):02}_00_{ew}{abs(lon):03}_00_DEM"
            )
            urls.append(
                f"https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif"
            )
    native = output / "glo30_native_window.tif"
    provenance_path = output / "reference.json"
    if native.exists() or provenance_path.exists():
        if not (native.exists() and provenance_path.exists()):
            raise ValueError(
                "Incomplete cached reference; inspect the output directory"
            )
        provenance = json.loads(provenance_path.read_text())
        if (
            provenance["bounds"] != list(bounds)
            or provenance["sha256"] != sha256(native)
            or not provenance.get("native_grid_preserved")
        ):
            raise ValueError(
                "Cached reference does not match this grid or its checksum"
            )
    else:
        with (
            rasterio.Env(
                GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                GDAL_HTTP_TIMEOUT="45",
                GDAL_HTTP_MAX_RETRY="2",
            ),
            ExitStack() as stack,
        ):
            datasets = []
            source_metadata = []
            for url in urls:
                print("Reading reference window:", url, flush=True)
                src = stack.enter_context(rasterio.open(url))
                if src.crs != rasterio.crs.CRS.from_epsg(4326):
                    raise ValueError(f"Unexpected source CRS: {src.crs}")
                datasets.append(src)
                source_metadata.append(
                    {
                        "url": url,
                        "crs": src.crs.to_string(),
                        "shape": list(src.shape),
                        "transform": list(src.transform),
                        "tags": src.tags(),
                        "nodata": src.nodata,
                    }
                )
            if any(src.res != datasets[0].res for src in datasets):
                raise ValueError("Reference tiles have differing resolutions")
            aligned_bounds = align_bounds(bounds, datasets[0].transform)
            raster, affine = merge(datasets, bounds=aligned_bounds, nodata=np.nan)
            profile = datasets[0].profile.copy()
            profile.update(
                height=raster.shape[1], width=raster.shape[2], transform=affine
            )
        write_grid(
            native,
            raster[0],
            profile,
            "Copernicus GLO-30 elevation (EGM2008 metres)",
            vertical_reference="EGM2008",
            source_release="Copernicus DEM 2021, AWS COG",
        )
        provenance = {
            "bounds": list(bounds),
            "sources": source_metadata,
            "fetched_utc": datetime.now(timezone.utc).isoformat(),
            "sha256": sha256(native),
            "vertical_reference": "EGM2008 (EPSG:3855)",
            "native_grid_preserved": True,
            "read_bounds": list(aligned_bounds),
            "source_release": "Copernicus DEM 2021, AWS COG",
            "documentation": "https://copernicus-dem-30m.s3.amazonaws.com/readme.html",
        }
        provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    with rasterio.open(native) as src:
        native_data = src.read(1, masked=True).filled(np.nan)
        native_transform = src.transform
        aligned = np.full(dem.shape, np.nan, dtype="float32")
        # Explicit 2-D horizontal CRS: no automatic vertical shift is requested.
        reproject(
            native_data,
            aligned,
            src_transform=src.transform,
            src_crs="EPSG:4326",
            src_nodata=np.nan,
            dst_transform=dem.transform,
            dst_crs=dem.crs,
            dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )
    return aligned, native_data, native_transform, provenance


def compare(dem_path, output, spots_path=None, distance_path=None):
    output.mkdir(parents=True, exist_ok=True)
    with rasterio.open(dem_path) as dem:
        soi = dem.read(1, masked=True).filled(np.nan)
        reference, native, native_transform, provenance = reference_window(dem, output)
        profile = dem.profile.copy()
        if not dem.crs.is_projected or dem.crs.linear_units != "metre":
            raise ValueError("Comparison expects a DEM projected in metres")
        dy, dx = abs(dem.transform.e), abs(dem.transform.a)
        gy, gx = np.gradient(reference, dy, dx)
        slope = np.degrees(np.arctan(np.hypot(gx, gy)))
        valid = np.isfinite(soi) & np.isfinite(reference)
        diff = np.where(valid, soi - reference, np.nan)
        interior = binary_erosion(valid, iterations=5)
        report = {
            "difference_sign": "SOI minus Copernicus GLO-30",
            "dem": str(dem_path),
            "dem_sha256": sha256(dem_path),
            "reference": provenance,
            "crs": dem.crs.to_string(),
            "shape": list(dem.shape),
            "transform": list(dem.transform),
            "resampling": "bilinear GLO-30 to the exact SOI pixel grid",
            "vertical_transformation_applied": False,
            "soi_vertical_reference": "Source VALUE unchanged; exact Indian vertical datum realization unverified",
            "common_area_km2": float(valid.sum() * dx * dy / 1e6),
            "soi_valid_pixels": int(np.isfinite(soi).sum()),
            "all_common_pixels": stats(diff),
            "interior_after_five_pixel_mask_erosion": stats(diff[interior]),
            "slope_groups": {},
            "limitations": [
                "Copernicus is a surface model including vegetation and buildings, not an independent bare-earth truth.",
                "Differences include interpolation, horizontal registration, source epoch and surface effects as well as vertical-reference differences.",
                "A sheet-wide offset does not identify the source vertical datum.",
                "Spatially correlated pixels are not independent accuracy observations.",
            ],
        }
        for lower, upper in [(0, 5), (5, 15), (15, 30), (30, 90)]:
            mask = interior & (slope >= lower) & (slope < upper)
            report["slope_groups"][f"{lower}-{upper}_degrees"] = stats(diff[mask])
        if distance_path:
            with rasterio.open(distance_path) as src:
                if (
                    src.shape != dem.shape
                    or src.transform != dem.transform
                    or src.crs != dem.crs
                ):
                    raise ValueError(
                        "Constraint-distance grid differs from the DEM grid"
                    )
                distance = src.read(1, masked=True).filled(np.nan)
            report["within_30m_of_elevation_constraint"] = stats(
                diff[interior & (distance <= 30)]
            )
            report["distance_grid_sha256"] = sha256(distance_path)
        point_rows = []
        if spots_path:
            features = [
                json.loads(line.lstrip("\x1e"))
                for line in spots_path.read_text().splitlines()
                if line.strip()
            ]
            lons = np.array([f["geometry"]["coordinates"][0] for f in features])
            lats = np.array([f["geometry"]["coordinates"][1] for f in features])
            sampled = sample_bilinear(native, native_transform, lons, lats)
            sx, sy = Transformer.from_crs(4326, dem.crs, always_xy=True).transform(
                lons, lats
            )
            point_slope = sample_bilinear(slope, dem.transform, sx, sy)
            for feature, z, s in zip(features, sampled, point_slope, strict=True):
                p = feature["properties"]
                point_rows.append(
                    {
                        "source_id": p.get("source_id"),
                        "symbol_type": p["symbol_type"],
                        "longitude": feature["geometry"]["coordinates"][0],
                        "latitude": feature["geometry"]["coordinates"][1],
                        "soi_m": p["elevation_m"],
                        "glo30_m": float(z),
                        "soi_minus_glo30_m": float(p["elevation_m"] - z),
                        "glo30_slope_degrees": float(s),
                    }
                )
            if point_rows:
                with (output / "point_differences.csv").open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(point_rows[0]))
                    writer.writeheader()
                    writer.writerows(point_rows)
                report["source_points"] = {
                    "sampling": "bilinear at exact source point coordinates in the native GLO-30 window",
                    "sha256": sha256(spots_path),
                    "all": stats([r["soi_minus_glo30_m"] for r in point_rows]),
                    "circle": stats(
                        [
                            r["soi_minus_glo30_m"]
                            for r in point_rows
                            if r["symbol_type"] == "circle"
                        ]
                    ),
                    "triangle": stats(
                        [
                            r["soi_minus_glo30_m"]
                            for r in point_rows
                            if r["symbol_type"] == "triangle"
                        ]
                    ),
                    "slope_below_5_degrees": stats(
                        [
                            r["soi_minus_glo30_m"]
                            for r in point_rows
                            if r["glo30_slope_degrees"] < 5
                        ]
                    ),
                }
    write_grid(
        output / "glo30_aligned.tif",
        reference,
        profile,
        "Copernicus GLO-30 elevation, EGM2008 metres",
        vertical_reference="EGM2008",
    )
    write_grid(
        output / "soi_minus_glo30.tif",
        diff,
        profile,
        "SOI minus Copernicus GLO-30 (metres)",
        vertical_shift_applied="false",
    )
    report["script_sha256"] = sha256(__file__)
    (output / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )

    extent = [0, soi.shape[1] * dx / 1000, 0, soi.shape[0] * dy / 1000]
    limits = np.percentile(reference[valid], [0, 100])
    color_limit = max(10, float(np.percentile(np.abs(diff[valid]), 98)))
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), layout="constrained")
    for ax, data, title in zip(
        axes.flat[:2],
        [soi, np.where(valid, reference, np.nan)],
        ["SOI contour + point interpolation", "Copernicus GLO-30 (EGM2008)"],
        strict=True,
    ):
        im = ax.imshow(
            data, extent=extent, cmap="terrain", vmin=limits[0], vmax=limits[1]
        )
        ax.set_title(title)
        fig.colorbar(im, ax=ax, label="Elevation (m)", shrink=0.85)
    ax = axes[1, 0]
    im = ax.imshow(
        diff, extent=extent, cmap="RdBu_r", vmin=-color_limit, vmax=color_limit
    )
    ax.set_title("SOI − GLO-30; red = SOI higher")
    fig.colorbar(im, ax=ax, label="Difference (m)", extend="both", shrink=0.85)
    for ax in [axes[0, 0], axes[0, 1], axes[1, 0]]:
        ax.set_xlabel("Easting from grid west edge (km)")
        ax.set_ylabel("Northing from grid south edge (km)")
    ax = axes[1, 1]
    ax.hist(
        diff[valid],
        bins=np.linspace(-color_limit, color_limit, 101),
        density=True,
        alpha=0.7,
    )
    ax.axvline(0, color="black", linewidth=1)
    ax.axvline(
        report["all_common_pixels"]["median_m"], color="firebrick", label="Median"
    )
    ax.set_xlabel("SOI − GLO-30 (m); tails outside plot omitted")
    ax.set_ylabel("Density within plotted range")
    ax.set_title("Raw differences; no fitted vertical correction")
    ax.legend()
    s = report["all_common_pixels"]
    fig.suptitle(
        f"{dem_path.parent.name}: mean {s['mean_m']:+.2f} m, median {s['median_m']:+.2f} m, RMSE {s['rmse_m']:.2f} m"
    )
    fig.savefig(output / "comparison.png", dpi=160)
    plt.close(fig)
    print(
        json.dumps(
            {
                "pixels": report["all_common_pixels"],
                "points": report.get("source_points"),
            },
            indent=2,
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dem", type=Path)
    parser.add_argument("--spots", type=Path)
    parser.add_argument("--constraint-distance", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    compare(args.dem, args.output_dir, args.spots, args.constraint_distance)
