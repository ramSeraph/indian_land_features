# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "rasterio>=1.4", "scipy>=1.14", "pyproj>=3.6", "shapely>=2.1", "matplotlib>=3.9"]
# ///
"""Validate local reconciliation on aligned neighbouring pilot sheets.

This intentionally rejects different UTM grids; general reprojection and
country-wide scheduling are not implied by this calibration experiment.
"""

import argparse
import hashlib
import importlib.metadata
import json
import math
import sys
import time
from functools import reduce
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.features import geometry_mask
from rasterio.shutil import copy as raster_copy
from rasterio.transform import from_origin
from shapely import segmentize
from shapely.geometry import mapping
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from contours_to_dem import get_sheet, transform, write_raster
from seam_interpolate import reconcile_patch


def file_record(path):
    path = Path(path)
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def run(sheet_dirs, output, replace_m=2000, context_m=3600):
    started = time.monotonic()
    if replace_m <= 0 or context_m < 3600:
        raise ValueError(
            "Use a positive replacement width and at least 3,600 m of context"
        )
    if len(sheet_dirs) < 2 or len(set(sheet_dirs)) != len(sheet_dirs):
        raise ValueError("Provide at least two distinct neighbouring sheets")
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Choose a new output directory")
    metadata, sheets = [], []
    for folder in sheet_dirs:
        report = json.loads((folder / "report.json").read_text())
        sheets.append(get_sheet(ROOT / "index_50k.geojson", report["sheet"]))
        with rasterio.open(folder / "dem.tif") as src:
            metadata.append((report, src.transform, src.crs, src.shape))
    crs, resolution = metadata[0][2], metadata[0][1].a
    if any(
        m[2] != crs or m[1].a != resolution or m[1].e != -resolution for m in metadata
    ):
        raise ValueError("Pilot requires one aligned metric grid")
    project = Transformer.from_crs(4326, crs, always_xy=True).transform
    joins = [
        a.boundary.intersection(b.boundary)
        for i, a in enumerate(sheets)
        for b in sheets[i + 1 :]
    ]
    shared = unary_union([g for g in joins if g.length > 0])
    if shared.is_empty:
        raise ValueError("Sheets have no shared edge")
    shared = transform(project, segmentize(shared, 0.002))
    domain = transform(project, segmentize(unary_union(sheets), 0.002))
    core_polygon = shared.buffer(replace_m).intersection(domain)
    context_polygon = shared.buffer(replace_m + context_m).intersection(domain)
    minx, miny, maxx, maxy = context_polygon.bounds
    west, south = (
        math.floor(minx / resolution) * resolution,
        math.floor(miny / resolution) * resolution,
    )
    east, north = (
        math.ceil(maxx / resolution) * resolution,
        math.ceil(maxy / resolution) * resolution,
    )
    affine = from_origin(west, north, resolution, resolution)
    grid_shape = (
        round((north - south) / resolution),
        round((east - west) / resolution),
    )
    inside = geometry_mask([mapping(domain)], grid_shape, affine, invert=True)
    core = (
        geometry_mask([mapping(core_polygon)], grid_shape, affine, invert=True) & inside
    )
    known, prior = np.full(grid_shape, np.nan), np.full(grid_shape, np.nan)
    unknown = np.zeros(grid_shape, dtype=np.int32)
    unknown_groups, offset = [], 0
    placements = []
    for folder, (report, old_transform, _, old_shape), polygon in zip(
        sheet_dirs, metadata, sheets, strict=True
    ):
        col, row = ~affine * (old_transform.c, old_transform.f)
        if abs(col - round(col)) > 1e-6 or abs(row - round(row)) > 1e-6:
            raise ValueError("Input pixel origins are not aligned")
        col, row = round(col), round(row)
        r0, r1 = max(0, row), min(grid_shape[0], row + old_shape[0])
        c0, c1 = max(0, col), min(grid_shape[1], col + old_shape[1])
        dst = np.s_[r0:r1, c0:c1]
        source = np.s_[r0 - row : r1 - row, c0 - col : c1 - col]
        placements.append((folder, report, old_transform, old_shape, dst, source))
        inside = geometry_mask(
            [mapping(transform(project, segmentize(polygon, 0.002)))],
            old_shape,
            old_transform,
            invert=True,
        )
        with rasterio.open(folder / "dem.tif") as src:
            data = src.read(1, masked=True).filled(np.nan)[source]
        with np.load(folder / "interpolation_inputs.npz") as inputs:
            k = inputs["known_heights"][source].copy()
            full_unknown = inputs["unknown_ids"].copy()
        full_unknown[~inside] = 0
        u = full_unknown[source].copy()
        full_core = np.zeros(old_shape, dtype=bool)
        full_core[source] = core[dst]
        for group in np.unique(full_unknown[full_core & (full_unknown > 0)]):
            if np.any((full_unknown == group) & ~full_core):
                raise ValueError(
                    f"{report['sheet']}: complete unknown group {group} crosses the replacement boundary; enlarge the strip"
                )
        # A raster-only pilot cannot establish that unknown fragments on the
        # two sheets represent the same source contour. Refuse that case until
        # their source-vector identities can be reconciled.
        near_join = geometry_mask(
            [mapping(shared.buffer(resolution))], old_shape, old_transform, invert=True
        )
        if np.any((full_unknown > 0) & near_join):
            raise ValueError(
                "Unknown contours touch a shared edge; reconcile their vector identities first"
            )
        k[~inside[source]] = np.nan
        u[~inside[source]] = 0
        both = np.isfinite(known[dst]) & np.isfinite(k)
        if np.any(both & (np.abs(known[dst] - k) > 1e-7)):
            raise ValueError(
                "Neighbouring constraint rasters disagree in a shared pixel"
            )
        known[dst][np.isfinite(k)] = k[np.isfinite(k)]
        prior[dst][np.isfinite(data)] = data[np.isfinite(data)]
        for group in np.unique(u[u > 0]):
            offset += 1
            unknown[dst][u == group] = offset
            unknown_groups.append(
                {
                    "joint_id": offset,
                    "sheet": report["sheet"],
                    "original_id": int(group),
                }
            )
    inside = geometry_mask([mapping(domain)], grid_shape, affine, invert=True)
    core = (
        geometry_mask([mapping(core_polygon)], grid_shape, affine, invert=True) & inside
    )
    patched, report = reconcile_patch(
        known, unknown, prior, core, progress=lambda stage: print(stage, flush=True)
    )
    patched[~inside] = np.nan
    output.mkdir(parents=True)
    write_raster(
        output / "before_patch.tif",
        prior,
        affine,
        crs.to_epsg(),
        "Independent sheet DEMs",
    )
    write_raster(
        output / "after_patch.tif",
        patched,
        affine,
        crs.to_epsg(),
        "Local shared-edge interpolation",
    )
    write_raster(
        output / "patch_change.tif",
        patched - prior,
        affine,
        crs.to_epsg(),
        "Reconciled minus independent elevation",
    )
    difference = np.abs(patched - prior)
    finite = core & np.isfinite(difference)
    report.update(
        sheets=[m[0]["sheet"] for m in metadata],
        replacement_width_per_side_m=replace_m,
        context_beyond_replacement_m=context_m,
        unknown_id_mapping=unknown_groups,
        common_crs=crs.to_string(),
        working_shape=list(grid_shape),
        changes={
            "count": int(finite.sum()),
            "mean_abs_m": float(np.mean(difference[finite])),
            "p99_abs_m": float(np.percentile(difference[finite], 99)),
            "max_abs_m": float(np.max(difference[finite])),
        },
        four_sheet_corners_included=len(sheets) == 4,
    )
    report["four_sheet_corners_included"] = (
        len(sheets) == 4
        and not reduce(
            lambda a, b: a.intersection(b), [s.boundary for s in sheets]
        ).is_empty
    )
    report["inputs"] = [
        file_record(folder / name)
        for folder in sheet_dirs
        for name in ("dem.tif", "report.json", "interpolation_inputs.npz")
    ]
    report["implementation"] = [
        file_record(ROOT / name)
        for name in (
            "index_50k.geojson",
            "imagico_interpolate.py",
            "seam_interpolate.py",
            "calibration/reconcile_adjacent_pilot.py",
        )
    ]
    report["environment"] = {
        "python": sys.version,
        "gdal": rasterio.__gdal_version__,
        **{
            name: importlib.metadata.version(name)
            for name in ("numpy", "scipy", "rasterio", "pyproj", "shapely")
        },
    }
    report["sheet_checks"] = {}
    report["limitations"].extend(
        [
            "This pilot requires aligned grids in the same UTM zone; nationwide scheduling is not implemented.",
            "Unknown contours that meet the shared edge are rejected pending vector identity reconciliation.",
            "Real contradictions between neighbouring source heights are not repaired by interpolating a shared strip.",
        ]
    )
    for folder, original, old_transform, old_shape, dst, source in placements:
        with rasterio.open(folder / "dem.tif") as src:
            data = src.read(1, masked=True).filled(np.nan)
        old = data.copy()
        local_core = core[dst]
        data[source][local_core] = patched[dst][local_core]
        # Each output retains its own sheet mask, including supported new cells.
        polygon = get_sheet(ROOT / "index_50k.geojson", original["sheet"])
        sheet_inside = geometry_mask(
            [mapping(transform(project, segmentize(polygon, 0.002)))],
            old_shape,
            old_transform,
            invert=True,
        )
        data[~sheet_inside] = np.nan
        with np.load(folder / "interpolation_inputs.npz") as inputs:
            measured = np.isfinite(inputs["known_heights"]) & sheet_inside
            if not np.array_equal(
                data[measured].astype("float32"),
                inputs["known_heights"][measured].astype("float32"),
            ):
                raise ValueError("A source constraint changed in the corrected sheet")
        full_core = np.zeros(old_shape, dtype=bool)
        full_core[source] = local_core
        if not np.array_equal(data[~full_core], old[~full_core], equal_nan=True):
            raise ValueError("Pixels outside the edge strip changed")
        plain = output / f"{original['sheet']}.working.tif"
        write_raster(
            plain,
            data,
            old_transform,
            crs.to_epsg(),
            "Elevation after local shared-edge pilot",
        )
        raster_copy(
            plain,
            output / f"{original['sheet']}.tif",
            driver="COG",
            compress="DEFLATE",
            predictor="FLOATING_POINT",
            overview_resampling="AVERAGE",
        )
        with rasterio.open(output / f"{original['sheet']}.tif") as src:
            saved = src.read(1, masked=True).filled(np.nan)
            if src.tags(ns="IMAGE_STRUCTURE").get(
                "LAYOUT"
            ) != "COG" or not np.array_equal(
                saved, data.astype("float32"), equal_nan=True
            ):
                raise ValueError("Corrected COG changed pixels or lacks COG layout")
        report["sheet_checks"][original["sheet"]] = {
            "measured_constraint_pixels_preserved": int(measured.sum()),
            "outside_strip_unchanged": True,
            "cog_pixels_identical": True,
        }
        plain.unlink()
    report["elapsed_seconds"] = time.monotonic() - started
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["changes"]), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sheets", nargs="+", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--replace-m", type=float, default=2000)
    p.add_argument("--context-m", type=float, default=3600)
    a = p.parse_args()
    run(a.sheets, a.output_dir, a.replace_m, a.context_m)
