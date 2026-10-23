# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "rasterio>=1.4", "pyproj>=3.6", "shapely>=2.1", "matplotlib>=3.9", "pillow>=11"]
# ///
"""Make raw map crops and source-vector overlays at audited contour conflicts.

The scanned sheet is calibration evidence only. This script never changes
contour geometry, source elevations, or DEM constraints.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from PIL import Image, ImageDraw, ImageFont
from pyproj import Transformer
from rasterio.windows import Window
from shapely.geometry import shape


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("toposheet", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--context-m", type=float, default=600)
    args = parser.parse_args()
    report = json.loads(args.report.read_text())
    source_lines = {}
    with (args.report.parent / "involved_source_lines.geojsonl").open() as stream:
        for line in stream:
            feature = json.loads(line)
            props = feature["properties"]
            source_lines[(props["source_id"], props["part_index"])] = feature
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows, cards = [], []
    font = ImageFont.load_default(size=17)
    with rasterio.open(args.toposheet) as raster:
        to_raster = Transformer.from_crs(4326, raster.crs, always_xy=True)
        to_utm = Transformer.from_crs(4326, report["projected_epsg"], always_xy=True)
        utm_to_raster = Transformer.from_crs(
            report["projected_epsg"], raster.crs, always_xy=True
        )
        for row in report["locations"]:
            lon, lat = row["coordinates"][:2]
            x, y = to_raster.transform(lon, lat)
            ux, uy = to_utm.transform(lon, lat)
            dx, dy = utm_to_raster.transform(ux + args.context_m, uy)
            width = max(32, math.ceil(math.hypot(dx - x, dy - y) / abs(raster.res[0])))
            px, py = ~raster.transform * (x, y)
            window = Window(round(px - width / 2), round(py - width / 2), width, width)
            data = raster.read([1, 2, 3], window=window, boundless=True, fill_value=255)
            data[:, raster.dataset_mask(window=window, boundless=True) == 0] = 255
            pixels = np.moveaxis(data, 0, 2)
            crop_transform = raster.window_transform(window)
            centre = np.array([px - window.col_off, py - window.row_off])
            cid = row["conflict_id"]
            Image.fromarray(pixels).save(args.output_dir / f"{cid}-raw-native.png")
            fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)
            for ax in axes:
                # Explicit pixel-boundary extent matches the affine transform.
                ax.imshow(pixels, extent=(0, width, width, 0), interpolation="nearest")
                ax.set(xlim=(0, width), ylim=(width, 0), aspect="equal")
                ax.set_xticks([])
                ax.set_yticks([])
                for delta in [(-6, 0), (6, 0), (0, -6), (0, 6)]:
                    start = centre + np.array(delta)
                    end = centre + np.array(delta) * 1.7
                    ax.plot(
                        [start[0], end[0]], [start[1], end[1]], color="#0072ff", lw=1.5
                    )
            axes[0].set_title("Scanned toposheet — blue crosshair marks conflict")
            axes[1].set_title("Same crop with original source contours")
            contacts = sorted(
                row["source_contacts"], key=lambda r: (r["elevation_m"], r["source_id"])
            )
            for color, contact in zip(
                ["#c000c0", "#00a000", "#e06000", "#0040ff"], contacts
            ):
                feature = source_lines[(contact["source_id"], contact["part_index"])]
                xy = np.asarray(shape(feature["geometry"]).coords)
                rx, ry = to_raster.transform(xy[:, 0], xy[:, 1])
                cx, cy = ~crop_transform * (rx, ry)
                label = f"{contact['elevation_m']:g} m / {contact['source_id'].split('.')[-1]}"
                axes[1].plot(cx, cy, color=color, lw=1.2, alpha=0.8, label=label)
            axes[1].legend(loc="upper right", fontsize=9)
            fig.suptitle(
                f"{report['sheet']} · {cid} · {lon:.7f}°, {lat:.7f}° · {args.context_m:g} m context"
            )
            fig.savefig(args.output_dir / f"{cid}-comparison.png", dpi=130)
            plt.close(fig)
            raw = Image.fromarray(pixels).resize((384, 384), Image.Resampling.NEAREST)
            marked = raw.copy()
            draw = ImageDraw.Draw(marked)
            mx, my = centre * 384 / width
            for a, b, c, d in [
                (-25, 0, -13, 0),
                (13, 0, 25, 0),
                (0, -25, 0, -13),
                (0, 13, 0, 25),
            ]:
                draw.line((mx + a, my + b, mx + c, my + d), fill=(0, 80, 255), width=2)
            card = Image.new("RGB", (384, 430), "white")
            card.paste(marked, (0, 46))
            title = cid + "   " + " / ".join(f"{v:g} m" for v in row["elevations_m"])
            ImageDraw.Draw(card).text((8, 8), title, fill="black", font=font)
            cards.append(card)
            rows.append(
                {
                    "conflict_id": cid,
                    "lon": lon,
                    "lat": lat,
                    "elevations_m": ";".join(str(v) for v in row["elevations_m"]),
                    "source_ids": ";".join(c["source_id"] for c in contacts),
                    "crop": f"{cid}-comparison.png",
                    "review": "",
                    "notes": "",
                }
            )
        raster_info = {
            "crs": str(raster.crs),
            "pixel_size": raster.res,
            "shape": raster.shape,
        }
    for start in range(0, len(cards), 9):
        canvas = Image.new("RGB", (384 * 3, 430 * 3), "#dddddd")
        for number, card in enumerate(cards[start : start + 9]):
            canvas.paste(card, ((number % 3) * 384, (number // 3) * 430))
        canvas.save(args.output_dir / f"contact-{start // 9 + 1:02}.png")
    review = args.output_dir / "review.csv"
    if not review.exists():
        with review.open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    manifest = {
        "sheet": report["sheet"],
        "conflict_locations": len(rows),
        "context_width_m": args.context_m,
        "raster": raster_info,
        "inputs": [
            {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in (args.report, args.toposheet)
        ],
        "note": "Manual calibration evidence; no automatic interpretation or source edits. Native crops preserve raster pixels.",
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
