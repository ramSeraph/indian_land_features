# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "rasterio>=1.4", "scipy>=1.14", "pyproj>=3.6", "shapely>=2.1", "matplotlib>=3.9"]
# ///
"""Compare a local seam pilot against cached, aligned GLO-30 reference windows."""

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.features import geometry_mask
from shapely import segmentize
from shapely.geometry import mapping
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from compare_dem import sha256, stats
from contours_to_dem import get_sheet, transform


def read(path):
    with rasterio.open(path) as src:
        return (
            src.read(1, masked=True).filled(np.nan).astype(float),
            src.transform,
            src.crs,
        )


def evaluate(patch_dir, reference_dir):
    patch_dir, reference_dir = Path(patch_dir), Path(reference_dir)
    provenance = json.loads((patch_dir / "report.json").read_text())
    before, affine, crs = read(patch_dir / "before_patch.tif")
    after, after_transform, after_crs = read(patch_dir / "after_patch.tif")
    if after_transform != affine or after_crs != crs or after.shape != before.shape:
        raise ValueError("Patch rasters differ in grid")
    reference = np.full(before.shape, np.nan)
    owners = np.zeros(before.shape, dtype=int)
    project = Transformer.from_crs(4326, crs, always_xy=True).transform
    polygons, inputs = [], []
    for owner, sheet in enumerate(provenance["sheets"], 1):
        polygon = get_sheet(ROOT / "index_50k.geojson", sheet)
        polygons.append(polygon)
        inside = geometry_mask(
            [mapping(transform(project, segmentize(polygon, 0.002)))],
            before.shape,
            affine,
            invert=True,
        )
        if np.any(inside & (owners > 0)):
            raise ValueError("Sheet interiors overlap")
        owners[inside] = owner
        path = reference_dir / sheet / "glo30_aligned.tif"
        data, old_transform, old_crs = read(path)
        if old_crs != crs or old_transform.a != affine.a or old_transform.e != affine.e:
            raise ValueError("Reference grid differs from patch grid")
        col, row = ~affine * (old_transform.c, old_transform.f)
        if abs(col - round(col)) > 1e-6 or abs(row - round(row)) > 1e-6:
            raise ValueError("Reference pixels are not aligned")
        col, row = round(col), round(row)
        r0, r1 = max(0, row), min(before.shape[0], row + data.shape[0])
        c0, c1 = max(0, col), min(before.shape[1], col + data.shape[1])
        dst = np.s_[r0:r1, c0:c1]
        source = np.s_[r0 - row : r1 - row, c0 - col : c1 - col]
        reference[dst][inside[dst]] = data[source][inside[dst]]
        inputs.append({"path": str(path), "sha256": sha256(path)})
    joins = [
        a.boundary.intersection(b.boundary)
        for i, a in enumerate(polygons)
        for b in polygons[i + 1 :]
    ]
    shared = transform(
        project, segmentize(unary_union([g for g in joins if g.length > 0]), 0.002)
    )
    common = np.isfinite(before) & np.isfinite(after) & np.isfinite(reference)
    bands = {}
    last = np.zeros(before.shape, bool)
    for width in sorted({300, 1000, provenance["replacement_width_per_side_m"]}):
        area = geometry_mask(
            [mapping(shared.buffer(width))], before.shape, affine, invert=True
        )
        selected = area & ~last & common
        bands[str(width)] = {
            "before": stats((before - reference)[selected]),
            "after": stats((after - reference)[selected]),
            "change": stats((after - before)[selected]),
        }
        last = area
    # Compare adjacent pixel steps across sheet boundaries. GLO-30's step is
    # subtracted so a real terrain slope is not mistaken for a sheet seam.
    steps = {
        "before": [],
        "after": [],
        "glo30": [],
        "before_minus_glo30": [],
        "after_minus_glo30": [],
    }
    for a, b in [(np.s_[:-1, :], np.s_[1:, :]), (np.s_[:, :-1], np.s_[:, 1:])]:
        selected = (
            common[a]
            & common[b]
            & (owners[a] > 0)
            & (owners[b] > 0)
            & (owners[a] != owners[b])
        )
        old_step = (before[b] - before[a])[selected]
        new_step = (after[b] - after[a])[selected]
        glo_step = (reference[b] - reference[a])[selected]
        for key, values in zip(
            steps,
            [old_step, new_step, glo_step, old_step - glo_step, new_step - glo_step],
            strict=True,
        ):
            steps[key].append(values)
    report = {
        "sheets": provenance["sheets"],
        "distance_bands_upper_bound_m": bands,
        "adjacent_pixel_steps_across_sheet_boundaries": {
            key: stats(np.concatenate(values)) for key, values in steps.items()
        },
        "reference_inputs": inputs,
        "patch_report_sha256": sha256(patch_dir / "report.json"),
        "script_sha256": sha256(__file__),
        "interpretation": "GLO-30 is a comparison surface, not a height input. Step residuals include terrain/reference mismatch; they are not survey accuracy estimates.",
    }
    transition_before, transition_after = [], []
    for a, b in [(np.s_[:-1, :], np.s_[1:, :]), (np.s_[:, :-1], np.s_[:, 1:])]:
        selected = common[a] & common[b] & (last[a] != last[b])
        transition_before.append((before[b] - before[a])[selected])
        transition_after.append((after[b] - after[a])[selected])
    transition_before = np.concatenate(transition_before)
    transition_after = np.concatenate(transition_after)
    report["replacement_boundary_pixel_steps"] = {
        "before": stats(transition_before),
        "after": stats(transition_after),
        "change": stats(transition_after - transition_before),
    }
    if provenance["four_sheet_corners_included"]:
        from functools import reduce

        corner = reduce(lambda a, b: a.intersection(b), [g.boundary for g in polygons])
        corner = transform(project, corner)
        selected = (
            geometry_mask(
                [mapping(corner.buffer(300))], before.shape, affine, invert=True
            )
            & common
        )
        report["within_300m_of_four_sheet_corner"] = {
            "before": stats((before - reference)[selected]),
            "after": stats((after - reference)[selected]),
        }
    (patch_dir / "glo30_check.json").write_text(json.dumps(report, indent=2) + "\n")
    valid_heights = np.concatenate([before[common], after[common], reference[common]])
    low, high = np.percentile(valid_heights, [1, 99])
    diff_limit = max(
        1.0, float(np.percentile(np.abs((before - reference)[common]), 98))
    )
    change_limit = max(
        1.0, float(np.percentile(np.abs((after - before)[common & last]), 99))
    )
    panels = [
        (before, "Independent sheets", "terrain", low, high),
        (after, "Local overlap interpolation", "terrain", low, high),
        (reference, "Copernicus GLO-30", "terrain", low, high),
        (before - reference, "Before − GLO-30", "RdBu_r", -diff_limit, diff_limit),
        (after - reference, "After − GLO-30", "RdBu_r", -diff_limit, diff_limit),
        (
            after - before,
            "Change from overlap pass",
            "RdBu_r",
            -change_limit,
            change_limit,
        ),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(18, 8), layout="constrained")
    extent = (
        affine.c / 1000,
        (affine.c + before.shape[1] * affine.a) / 1000,
        (affine.f + before.shape[0] * affine.e) / 1000,
        affine.f / 1000,
    )
    for ax, (values, title, cmap, vmin, vmax) in zip(axes.flat, panels, strict=True):
        image = ax.imshow(
            np.where(owners > 0, values, np.nan),
            extent=extent,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        lines = list(shared.geoms) if hasattr(shared, "geoms") else [shared]
        for line in lines:
            xy = np.asarray(line.coords) / 1000
            ax.plot(xy[:, 0], xy[:, 1], color="magenta", linewidth=0.7)
        ax.set_title(title)
        ax.set_xlabel("UTM easting (km)")
        ax.set_ylabel("UTM northing (km)")
        fig.colorbar(image, ax=ax, shrink=0.8, label="m")
    fig.suptitle(" / ".join(provenance["sheets"]) + " — magenta: shared sheet boundary")
    fig.savefig(patch_dir / "overlap_comparison.png", dpi=170)
    plt.close(fig)
    print(
        json.dumps(
            {
                "bands": bands,
                "steps": report["adjacent_pixel_steps_across_sheet_boundaries"],
            },
            indent=2,
        )
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("patch_dir", type=Path)
    parser.add_argument(
        "--reference-dir", type=Path, default=Path("data/comparisons/production-pilot")
    )
    args = parser.parse_args()
    evaluate(args.patch_dir, args.reference_dir)
