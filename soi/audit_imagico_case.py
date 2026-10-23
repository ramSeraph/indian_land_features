# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "rasterio>=1.4", "scipy>=1.14", "shapely>=2.1"]
# ///
# SPDX-License-Identifier: GPL-3.0-or-later
"""Calibrate unknown-height interpolation on a vector neighbourhood.

Explicitly selected source IDs are a calibration experiment, not an automatic
error detector. No toposheet or reference DEM is read. Inputs stay unchanged.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from pyproj import Transformer
from rasterio.transform import from_origin
from shapely.geometry import box, shape
from shapely.ops import unary_union

import contours_to_dem as dem
from imagico_dem import build_grid, sample_grid


def run_case(
    classified_dir,
    sheet_id,
    index,
    source_ids,
    output_dir,
    context_m=1200,
    resolution=5,
):
    if not all(math.isfinite(v) and v > 0 for v in (context_m, resolution)):
        raise ValueError("Context and resolution must be positive and finite")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("Output directory must be empty")
    sheet = dem.get_sheet(index, sheet_id)
    zone = min(60, int((sheet.centroid.x + 180) // 6) + 1)
    epsg = (32600 if sheet.centroid.y >= 0 else 32700) + zone
    project = Transformer.from_crs(4326, epsg, always_xy=True).transform
    unproject = Transformer.from_crs(epsg, 4326, always_xy=True).transform

    def load(name):
        return [
            (dem.transform(project, shape(f["geometry"])), f["properties"])
            for f in dem.read_features(classified_dir / name)
        ]

    contours, spots = load("contours.geojsonl"), load("spots.geojsonl")
    groups = {source_id: i + 1 for i, source_id in enumerate(dict.fromkeys(source_ids))}
    selected = [(g, p) for g, p in contours if p["source_id"] in groups]
    missing = set(groups) - {p["source_id"] for _, p in selected}
    if missing:
        raise ValueError(f"Unknown source IDs: {sorted(missing)}")
    centre = unary_union([g for g, _ in selected]).centroid
    west = math.floor((centre.x - context_m / 2) / resolution) * resolution
    north = math.ceil((centre.y + context_m / 2) / resolution) * resolution
    size = math.ceil(context_m / resolution)
    if size < 3 or size * size > 1_000_000:
        raise ValueError("Use a calibration window between 3x3 and 1 million pixels")
    area = box(west, north - size * resolution, west + size * resolution, north)
    grid_transform = from_origin(west, north, resolution, resolution)
    known, unknown = [], []
    for line, props in contours:
        for part in dem.parts(line.intersection(area), "LineString"):
            if part.length == 0:
                continue
            if props["source_id"] in groups:
                unknown.append(
                    (
                        part,
                        dict(
                            props,
                            VALUE=None,
                            original_elevation_m=props["VALUE"],
                            height_status="unassigned",
                            unknown_height_id=groups[props["source_id"]],
                            selection_method="explicit_calibration_selection",
                        ),
                    )
                )
            else:
                known.append((part, props))
    points = [(g, p) for g, p in spots if area.covers(g)]
    output_dir.mkdir(parents=True, exist_ok=True)
    xy, z = dem.constraint_samples(known, resolution / 2)
    if points:
        xy = np.concatenate([xy, [[g.x, g.y] for g, _ in points]])
        z = np.concatenate([z, [p["VALUE"] for _, p in points]])
    dem.unique_samples(
        xy,
        z,
        conflict_path=output_dir / "conflicting_constraints.geojsonl",
        unproject=unproject,
    )
    for name, records in [
        ("known_contours", known),
        ("unknown_height_contours", unknown),
        ("spots", points),
    ]:
        dem.write_vectors(output_dir / f"{name}.geojsonl", records, unproject)
    validation, _, _ = build_grid(
        known,
        [],
        unknown,
        (size, size),
        grid_transform,
        progress=lambda stage: print(f"Validation: {stage}", flush=True),
    )
    held_out = (
        dem.error_stats(
            sample_grid(
                validation.elevation, grid_transform, [[g.x, g.y] for g, _ in points]
            )
            - [p["VALUE"] for _, p in points]
        )
        if points
        else {"count": 0, "unavailable": 0}
    )
    dem.write_raster(
        output_dir / "dem_contours_only.tif",
        validation.elevation,
        grid_transform,
        epsg,
        "Validation only: point heights held out",
        method="imagico_python_holdout",
    )
    validation_report = dict(
        validation.report,
        unknown_height_estimates=validation.unknown_heights,
        point_errors=held_out,
    )
    if points:
        del validation
        result, known_grid, unknown_grid = build_grid(
            known,
            points,
            unknown,
            (size, size),
            grid_transform,
            progress=lambda stage: print(f"Final, all points: {stage}", flush=True),
        )
    else:
        result = validation
        # There are no point constraints to add; both passes have identical inputs.
        from imagico_dem import rasterize_constraints

        known_grid, unknown_grid, _ = rasterize_constraints(
            known, points, unknown, (size, size), grid_transform
        )
    dem.write_raster(
        output_dir / "dem.tif",
        result.elevation,
        grid_transform,
        epsg,
        "Calibration DEM with all local height points",
        method="imagico_python",
    )
    np.savez_compressed(
        output_dir / "interpolation_inputs.npz",
        known_heights=known_grid,
        unknown_ids=unknown_grid,
    )
    estimates = {p["unknown_id"]: p for p in result.unknown_heights}
    dem.write_vectors(
        output_dir / "estimated_contours.geojsonl",
        [
            (
                g,
                dict(
                    p,
                    VALUE=estimates[p["unknown_height_id"]]["estimated_elevation_m"],
                    height_status="estimated",
                    height_estimation=estimates[p["unknown_height_id"]],
                ),
            )
            for g, p in unknown
        ],
        unproject,
    )
    report = {
        "sheet": sheet_id,
        "method": "imagico_python",
        "selection_method": "explicit_calibration_selection",
        "automatically_selected": False,
        "unknown_source_ids": list(groups),
        "original_unknown_elevations_m": sorted(
            {p["original_elevation_m"] for _, p in unknown}
        ),
        "unknown_length_m": sum(g.length for g, _ in unknown),
        "context_m": context_m,
        "resolution_m": resolution,
        "shape": [size, size],
        "crs": f"EPSG:{epsg}",
        "point_heights_used": len(points),
        "source_files_changed": False,
        "validation": validation_report,
        "final": dict(result.report, unknown_height_estimates=result.unknown_heights),
        "inputs": [
            {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for path in [
                classified_dir / "contours.geojsonl",
                classified_dir / "spots.geojsonl",
                index,
            ]
        ],
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {k: report[k] for k in ("sheet", "resolution_m", "point_heights_used")}
            | {"estimates": result.unknown_heights},
            indent=2,
        )
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("classified_dir", type=Path)
    parser.add_argument("--sheet", required=True)
    parser.add_argument(
        "--index", type=Path, default=Path(__file__).with_name("index_50k.geojson")
    )
    parser.add_argument("--unknown-source", action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--context-m", type=float, default=1200)
    parser.add_argument("--resolution", type=float, default=5)
    args = parser.parse_args()
    try:
        run_case(
            args.classified_dir,
            args.sheet,
            args.index,
            args.unknown_source,
            args.output_dir,
            args.context_m,
            args.resolution,
        )
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
