# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "rasterio>=1.4", "scipy>=1.14", "pyproj>=3.6", "matplotlib>=3.9", "shapely>=2.1"]
# ///
"""Inspect native GLO-30 elevations around a contour conflict, without a SOI DEM."""

import argparse
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.coords import BoundingBox
from rasterio.transform import from_origin
from shapely.geometry import shape

from compare_dem import reference_window, sample_bilinear


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("contours", type=Path)
    parser.add_argument("--case", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--context-m", type=float, default=2000)
    args = parser.parse_args()
    if not math.isfinite(args.context_m) or args.context_m <= 0:
        parser.error("--context-m must be finite and positive")
    report = json.loads(args.report.read_text())
    case = next(r for r in report["locations"] if r["conflict_id"] == args.case)
    lon, lat = case["coordinates"][:2]
    crs = rasterio.crs.CRS.from_epsg(report["projected_epsg"])
    project = Transformer.from_crs(4326, crs, always_xy=True)
    unproject = Transformer.from_crs(crs, 4326, always_xy=True)
    cx, cy = project.transform(lon, lat)
    n = math.ceil(args.context_m / 30)
    half = n * 30 / 2
    grid = SimpleNamespace(
        crs=crs,
        shape=(n, n),
        bounds=BoundingBox(cx - half, cy - half, cx + half, cy + half),
        transform=from_origin(cx - half, cy + half, 30, 30),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    _, native, affine, _provenance = reference_window(grid, args.output_dir)
    rr, cc = np.indices(native.shape)
    xx, yy = affine * (cc + 0.5, rr + 0.5)
    ex, ny = project.transform(xx, yy)
    distance = np.hypot(ex - cx, ny - cy)
    point_height = float(sample_bilinear(native, affine, [lon], [lat])[0])
    column, row = ~affine * (lon, lat)
    pixel_height = float(native[math.floor(row), math.floor(column)])
    neighbourhoods = {}
    for radius in [30, 50, 100, 250, 500, 1000]:
        values = native[(distance <= radius) & np.isfinite(native)]
        neighbourhoods[str(radius)] = {
            "radius_m": radius,
            "native_pixel_centres": int(values.size),
            "min_m": float(values.min()) if values.size else None,
            "mean_m": float(values.mean()) if values.size else None,
            "max_m": float(values.max()) if values.size else None,
        }
    target_ids = {r["source_id"] for r in case["source_contacts"]}
    targets = []
    with args.contours.open() as stream:
        for line in stream:
            feature = json.loads(line)
            props = feature["properties"]
            if props.get("source_id", props.get("id")) in target_ids:
                coords = np.asarray(shape(feature["geometry"]).coords)
                z = sample_bilinear(native, affine, coords[:, 0], coords[:, 1])
                tx, ty = project.transform(coords[:, 0], coords[:, 1])
                targets.append((np.asarray(tx) - cx, np.asarray(ty) - cy, props, z))
    summary = {
        "sheet": report["sheet"],
        "case": args.case,
        "coordinates": [lon, lat],
        "reference": "Standard Copernicus GLO-30 AWS COG",
        "glo30_at_conflict_bilinear_m": point_height,
        "glo30_containing_native_pixel_m": pixel_height,
        "native_neighbourhoods": neighbourhoods,
        "source_fragment_samples": [
            {
                "source_id": p.get("source_id", p.get("id")),
                "soi_height_m": p["VALUE"],
                "glo30_at_source_vertices_m": z.tolist(),
            }
            for _, _, p, z in targets
        ],
        "comparison_scope": "Local terrain shape and evidence for a 1240 m feature, not a vertical-datum estimate.",
        "native_grid_preserved": True,
        "reference_provenance": "reference.json",
        "limitation": "GLO-30 smooths small terrain features. It cannot exclude a feature substantially narrower than its roughly 30 m sampling, and slope-dependent differences are expected.",
        "input_report": str(args.report),
        "input_contours": str(args.contours),
        "contours_sha256": hashlib.sha256(args.contours.read_bytes()).hexdigest(),
    }
    (args.output_dir / "report.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    with rasterio.open(args.output_dir / "glo30_native_window.tif") as src:
        rx, ry = src.res
    xstep = math.hypot(*np.subtract(project.transform(lon + rx, lat), (cx, cy)))
    ystep = math.hypot(*np.subtract(project.transform(lon, lat + ry), (cx, cy)))

    fig = plt.figure(figsize=(15, 11), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.65])
    ax0, ax1, ax2 = (
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, :]),
    )
    for ax, limit, title in [
        (ax0, 500, "1 km neighbourhood"),
        (ax1, 120, "240 m close-up — native pixels"),
    ]:
        mask = (np.abs(ex - cx) <= limit) & (np.abs(ny - cy) <= limit)
        mesh = ax.pcolormesh(
            ex - cx,
            ny - cy,
            native,
            shading="nearest",
            cmap="terrain",
            vmin=float(native[mask].min()),
            vmax=float(native[mask].max()),
        )
        fig.colorbar(mesh, ax=ax, label="GLO-30 height (m)", shrink=0.9)
        ax.set(
            xlim=(-limit, limit),
            ylim=(-limit, limit),
            aspect="equal",
            xlabel="East of conflict (m)",
            ylabel="North of conflict (m)",
            title=title,
        )
        for tx, ty, props, _ in targets:
            height = float(props["VALUE"])
            ax.plot(
                tx,
                ty,
                color="#d00000" if height == max(case["elevations_m"]) else "#111111",
                lw=2.2,
                label=f"SOI {height:g} m",
            )
        ax.plot(0, 0, "+", color="white", ms=14, mew=2.5)
        ax.plot(0, 0, "+", color="black", ms=12, mew=1)
        ax.legend(loc="upper left", fontsize=9)
        if limit == 120:
            for r, c in zip(*np.where(mask)):
                ax.text(
                    ex[r, c] - cx,
                    ny[r, c] - cy,
                    f"{native[r, c]:.0f}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="black",
                    bbox={"facecolor": "white", "alpha": 0.5, "edgecolor": "none", "pad": 0.2},
                )
    offsets = np.arange(-500, 501, 5)
    for orientation, xs, ys in [
        ("West → east", cx + offsets, np.full_like(offsets, cy, dtype=float)),
        ("South → north", np.full_like(offsets, cx, dtype=float), cy + offsets),
    ]:
        px, py = unproject.transform(xs, ys)
        values = sample_bilinear(native, affine, px, py)
        ax2.plot(offsets, values, label=orientation, lw=2)
    for height in case["elevations_m"]:
        ax2.axhline(
            height,
            color="#c00000" if height == max(case["elevations_m"]) else "#333333",
            ls="--",
            label=f"SOI {height:g} m",
        )
    ax2.axvline(0, color="#888888", lw=1)
    ax2.set(
        xlabel="Distance from conflict (m)",
        ylabel="Height (m)",
        title="GLO-30 profiles through C005; zero marks the conflict",
    )
    ax2.grid(alpha=0.2)
    ax2.legend(loc="upper left", ncols=4)
    fig.suptitle(
        f"{report['sheet']} · {args.case}/C006 · Copernicus GLO-30\n"
        f"At C005: {point_height:.1f} m · native spacing here ≈ {xstep:.1f} × {ystep:.1f} m",
        fontsize=16,
    )
    fig.savefig(args.output_dir / "terrain-check.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
