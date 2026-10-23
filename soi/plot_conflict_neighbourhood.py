# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "rasterio>=1.4", "pyproj>=3.6", "shapely>=2.1", "matplotlib>=3.9", "pillow>=11"]
# ///
"""Plot nearby source contours by recorded height for a conflict review."""

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import BoundaryNorm, ListedColormap
from PIL import Image
from pyproj import Transformer
from rasterio.windows import Window
from shapely.geometry import LineString, box, shape


def parts(geometry):
    if geometry.geom_type == "LineString":
        yield geometry
    elif hasattr(geometry, "geoms"):
        for component in geometry.geoms:
            yield from parts(component)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("contours", type=Path)
    parser.add_argument("toposheet", type=Path)
    parser.add_argument("--case", required=True)
    parser.add_argument("--context-m", type=float, default=900)
    parser.add_argument("--colour-levels", type=float, nargs="+",
                        help="Use a shared set of height colours across multiple views")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.context_m) or args.context_m <= 0:
        parser.error("--context-m must be finite and positive")
    report = json.loads(args.report.read_text())
    case = next(r for r in report["locations"] if r["conflict_id"] == args.case)
    lon, lat = case["coordinates"][:2]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{args.case}-height-context-{args.context_m:g}m"
    visible = []
    with rasterio.open(args.toposheet) as raster:
        to_raster = Transformer.from_crs(4326, raster.crs, always_xy=True)
        to_utm = Transformer.from_crs(4326, report["projected_epsg"], always_xy=True)
        utm_to_raster = Transformer.from_crs(
            report["projected_epsg"], raster.crs, always_xy=True
        )
        x, y = to_raster.transform(lon, lat)
        ux, uy = to_utm.transform(lon, lat)
        dx, dy = utm_to_raster.transform(ux + args.context_m, uy)
        width = math.ceil(math.hypot(dx - x, dy - y) / abs(raster.res[0]))
        px, py = ~raster.transform * (x, y)
        window = Window(round(px - width / 2), round(py - width / 2), width, width)
        data = raster.read([1, 2, 3], window=window, boundless=True, fill_value=255)
        data[:, raster.dataset_mask(window=window, boundless=True) == 0] = 255
        pixels = np.moveaxis(data, 0, 2)
        crop_transform = raster.window_transform(window)
        crop_bounds = box(0, 0, width, width)
        with args.contours.open() as stream:
            for text in stream:
                if not text.strip():
                    continue
                feature = json.loads(text.strip().lstrip("\x1e"))
                props = feature["properties"]
                height = float(props["VALUE"])
                for component in parts(shape(feature["geometry"])):
                    coords = np.asarray(component.coords)
                    rx, ry = to_raster.transform(coords[:, 0], coords[:, 1])
                    cx, cy = ~crop_transform * (rx, ry)
                    clipped = LineString(np.column_stack([cx, cy])).intersection(crop_bounds)
                    for line in parts(clipped):
                        if line.length > 0:
                            visible.append((height, line, props.get("source_id", props.get("id"))))

        locations = []
        for row in report["locations"]:
            rx, ry = to_raster.transform(*row["coordinates"][:2])
            cx, cy = ~crop_transform * (rx, ry)
            if 0 < cx < width and 0 < cy < width:
                locations.append((row["conflict_id"], cx, cy))

    visible_heights = sorted({height for height, _, _ in visible})
    if not visible_heights:
        raise ValueError("No contours in the requested crop")
    heights = sorted(set(args.colour_levels or visible_heights))
    if not set(visible_heights).issubset(heights) or not all(math.isfinite(h) for h in heights):
        parser.error("--colour-levels must be finite and include every visible height")
    colours = matplotlib.colormaps["turbo"](np.linspace(0.04, 0.96, len(heights)))
    by_height = dict(zip(heights, colours))
    cmap = ListedColormap(colours)
    norm = BoundaryNorm(np.arange(len(heights) + 1) - 0.5, len(heights))
    fig, axes = plt.subplots(1, 2, figsize=(16, 8.7), layout="constrained")
    axes[0].imshow(pixels, extent=(0, width, width, 0), interpolation="nearest")
    axes[0].imshow(np.ones_like(pixels) * 255, extent=(0, width, width, 0), alpha=0.37)
    axes[0].set_title("Source contours over the scanned toposheet", fontsize=14)
    axes[1].set_title("Same contours without the raster", fontsize=14)
    for ax in axes:
        for height, line, _ in visible:
            xy = np.asarray(line.coords)
            ax.plot(xy[:, 0], xy[:, 1], color=by_height[height], lw=1.8, solid_capstyle="round")
        for cid, cx, cy in locations:
            ax.plot(cx, cy, "o", ms=9 if cid == args.case else 6,
                    mfc="none", mec="black", mew=1.3, zorder=8)
            ax.annotate(cid, (cx, cy), xytext=(10, -18 if cid == args.case else 12),
                        textcoords="offset points", fontsize=10, weight="bold", zorder=10,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=0.88, pad=1.5))
        ax.set(xlim=(0, width), ylim=(width, 0), aspect="equal")
        ax.set_xticks([])
        ax.set_yticks([])

    # Label long visible pieces on the uncluttered panel, avoiding other labels
    # and the conflict markers. All values come directly from source VALUE.
    occupied = [(cx, cy) for _, cx, cy in locations]
    label_counts = {}
    for height, line, _ in sorted(visible, key=lambda record: -record[1].length):
        if line.length < width * 0.10 or label_counts.get(height, 0) >= 3:
            continue
        for fraction in (0.50, 0.30, 0.70, 0.15, 0.85):
            point = line.interpolate(fraction, normalized=True)
            tx, ty = point.x, point.y
            if not (width * 0.05 < tx < width * 0.95 and width * 0.04 < ty < width * 0.96):
                continue
            if any(abs(tx - ox) < width * 0.075 and abs(ty - oy) < width * 0.042 for ox, oy in occupied):
                continue
            axes[1].text(tx, ty, f"{height:g}", fontsize=9, ha="center", va="center",
                         color="#111111", zorder=6,
                         path_effects=[pe.withStroke(linewidth=3.8, foreground="white")])
            occupied.append((tx, ty))
            label_counts[height] = label_counts.get(height, 0) + 1
            break
    # A very short erroneous fragment can disappear under the conflict markers.
    # Give it an explicit height callout rather than relying on its colour alone.
    case_source_ids = {row["source_id"] for row in case["source_contacts"]}
    for height, line, source_id in visible:
        if source_id in case_source_ids and line.length < width * 0.10:
            point = line.interpolate(0.5, normalized=True)
            axes[1].annotate(
                f"{height:g} m", (point.x, point.y), xytext=(45, -45),
                textcoords="offset points", color=by_height[height], fontsize=11,
                weight="bold", zorder=11,
                arrowprops=dict(arrowstyle="->", color=by_height[height], lw=1.3),
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.94, pad=2),
            )
    bar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap),
                       ax=axes, ticks=np.arange(len(heights)), fraction=0.035, pad=0.015,
                       shrink=0.88)
    bar.ax.set_yticklabels([f"{height:g}" for height in heights])
    bar.set_label("Recorded contour height (m)", fontsize=12)
    fig.suptitle(f"{report['sheet']} · {args.case} · {args.context_m:g} m neighbourhood\n"
                 f"{lon:.7f}° E, {lat:.7f}° N · circles mark height conflicts", fontsize=16)
    fig.savefig(args.output_dir / f"{stem}.png", dpi=160)
    fig.savefig(args.output_dir / f"{stem}.svg")
    plt.close(fig)
    Image.fromarray(pixels).save(args.output_dir / f"{stem}-raw.png")
    manifest = {
        "sheet": report["sheet"], "case": args.case, "coordinates": [lon, lat],
        "context_width_m": args.context_m, "height_levels_m": visible_heights,
        "colour_levels_m": heights,
        "visible_parts": len(visible), "nearby_conflict_ids": [r[0] for r in locations],
        "inputs": {"report": str(args.report), "contours": str(args.contours),
                   "toposheet": str(args.toposheet)},
        "source_ids": sorted({source_id for _, _, source_id in visible if source_id}),
        "note": "Colours use recorded VALUE. Clipped only for display; no source edits. Raster faded in overlay; raw native crop saved separately.",
    }
    (args.output_dir / f"{stem}.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"image": str(args.output_dir / f"{stem}.png"), **manifest}, indent=2))


if __name__ == "__main__":
    main()
