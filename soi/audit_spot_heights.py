#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "rasterio>=1.4", "scipy>=1.14", "shapely>=2.1", "pillow>=11"]
# ///
"""Make reproducible map crops for visual review of the spot-symbol heuristic."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import rasterio
from PIL import Image, ImageDraw, ImageFont
from pyproj import Transformer
from rasterio.windows import Window
from shapely import normalize
from shapely.geometry import shape

from contours_to_dem import (
    geometry_hash,
    get_sheet,
    read_features,
    reconstruct_spot_symbols,
    roundness,
    spot_ring,
    spot_triangle,
    transform,
)


def audit(
    source,
    toposheet,
    index,
    sheet_id,
    output_dir,
    context_m=600,
    triangle_side=70.0,
    only_ids=(),
):
    output = Path(output_dir)
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    sheet = get_sheet(index, sheet_id)
    utm = (32600 if sheet.centroid.y >= 0 else 32700) + min(
        60, int((sheet.centroid.x + 180) // 6) + 1
    )
    project = Transformer.from_crs(4326, utm, always_xy=True).transform
    unproject = Transformer.from_crs(utm, 4326, always_xy=True).transform
    candidates, controls, seen = [], [], set()
    reconstructed = reconstruct_spot_symbols(
        read_features(source), project, triangle_side_m=triangle_side
    )
    requested_ids, matched_ids = set(only_ids), set()
    for feature in read_features(source):
        source_id = feature["properties"].get("id", "")
        if requested_ids and source_id not in requested_ids:
            continue
        matched_ids.add(source_id)
        geom = shape(feature["geometry"])
        components = geom.geoms if geom.geom_type == "MultiLineString" else [geom]
        for part_index, component in enumerate(components):
            value = float(feature["properties"]["VALUE"])
            key = normalize(component).wkb, value
            source_fragments = None
            if key in reconstructed:
                component, source_fragments = reconstructed[key]
                key = normalize(component).wkb, value
            if key in seen:
                continue
            seen.add(key)
            projected = transform(project, component)
            ring = spot_ring(projected)
            triangle = spot_triangle(projected, triangle_side) if ring is None else None
            centre = (
                ring[0] if ring else triangle[0] if triangle else projected.centroid
            )
            lon, lat = unproject(centre.x, centre.y)
            if not sheet.covers(transform(unproject, centre)):
                continue
            row = {
                "source_id": feature["properties"].get("id", ""),
                "part_index": part_index,
                "value_m": value,
                "lon": lon,
                "lat": lat,
                "length_m": projected.length,
                "vertices": len(component.coords),
                "closed": projected.is_closed,
                "radius_m": ring[1] if ring else None,
                "side_min_m": triangle[1]["symbol_side_min_m"] if triangle else None,
                "side_max_m": triangle[1]["symbol_side_max_m"] if triangle else None,
                "edge_error_m": triangle[1]["symbol_edge_error_m"]
                if triangle
                else None,
                "roundness": roundness(projected),
                "geometry_sha256": geometry_hash(component),
                "source_fragment_count": len(source_fragments)
                if source_fragments
                else 1,
                "group": "circle" if ring else "triangle" if triangle else "control",
            }
            if ring or triangle:
                candidates.append(row)
            elif projected.is_closed or projected.length < 250:
                controls.append(row)
    candidates.sort(
        key=lambda row: (
            row["group"] != "circle",
            (row["radius_m"] or 0) > 10,
            -row["lat"],
            row["lon"],
        )
    )
    if requested_ids - matched_ids:
        raise ValueError(f"Source IDs not found: {sorted(requested_ids - matched_ids)}")
    # Prefer actual small closed contours over sub-metre fragments at tile edges.
    controls.sort(
        key=lambda row: (not row["closed"], row["length_m"] < 1, row["length_m"])
    )
    # The smallest rejected geometries are the useful near-miss comparison set.
    rows = candidates + controls[:24]
    if not rows:
        raise ValueError(
            "No symbol candidates or small comparison contours in this selection"
        )
    font = ImageFont.load_default(size=16)
    cards = []
    with rasterio.open(toposheet) as raster:
        to_raster = Transformer.from_crs(4326, raster.crs, always_xy=True).transform
        from_utm = Transformer.from_crs(utm, raster.crs, always_xy=True).transform
        for number, row in enumerate(rows, 1):
            row["audit_id"] = f"{number:03d}"
            x, y = to_raster(row["lon"], row["lat"])
            px, py = (~raster.transform) * (x, y)
            ux, uy = project(row["lon"], row["lat"])
            rx, ry = from_utm(ux + context_m, uy)
            width = max(32, round(np.hypot(rx - x, ry - y) / abs(raster.res[0])))
            window = Window(round(px - width / 2), round(py - width / 2), width, width)
            data = raster.read([1, 2, 3], window=window, boundless=True, fill_value=255)
            valid = raster.dataset_mask(window=window, boundless=True)
            data[:, valid == 0] = 255
            raw = Image.fromarray(np.moveaxis(data, 0, 2)).resize((384, 384))
            raw.save(output / f"{number:03d}-raw.png")
            marked = raw.copy()
            draw = ImageDraw.Draw(marked)
            cx, cy = (
                (px - window.col_off) * 384 / width,
                (py - window.row_off) * 384 / width,
            )
            # Leave the symbol itself unobscured.
            for a, b, c, d in [
                (-22, 0, -12, 0),
                (12, 0, 22, 0),
                (0, -22, 0, -12),
                (0, 12, 0, 22),
            ]:
                draw.line((cx + a, cy + b, cx + c, cy + d), fill=(0, 110, 255), width=2)
            marked.save(output / f"{number:03d}-marked.png")
            card = Image.new("RGB", (384, 434), "white")
            card.paste(marked, (0, 50))
            draw = ImageDraw.Draw(card)
            detail = (
                f"r={row['radius_m']:.2f}m"
                if row["radius_m"]
                else f"sides={row['side_min_m']:.1f}-{row['side_max_m']:.1f}m"
                if row["side_min_m"]
                else f"{row['vertices']} vertices; {row['length_m']:.0f}m long"
            )
            draw.text(
                (8, 3),
                f"{number:03d}  {row['group']}  VALUE={row['value_m']:.3f}m",
                font=font,
                fill="black",
            )
            draw.text((8, 25), detail, font=font, fill="black")
            cards.append(card)
        raster_info = {
            "crs": str(raster.crs),
            "pixel_size": raster.res,
            "width": raster.width,
            "height": raster.height,
        }
    for start in range(0, len(cards), 12):
        contact = Image.new("RGB", (384 * 4, 434 * 3), "#dddddd")
        for position, card in enumerate(cards[start : start + 12]):
            contact.paste(card, ((position % 4) * 384, (position // 4) * 434))
        contact.save(output / f"contact-{start // 12 + 1:02d}.jpg", quality=95)
    with (output / "candidates.csv").open("w") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]) + ["review", "printed_height_m", "notes"]
        )
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "sheet": sheet_id,
        "context_width_m": context_m,
        "raster": raster_info,
        "files": [
            {
                "path": str(path),
                "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            }
            for path in (source, toposheet)
        ],
        "circle_candidates": sum(row["group"] == "circle" for row in candidates),
        "triangle_candidates": sum(row["group"] == "triangle" for row in candidates),
        "triangle_nominal_side_m": triangle_side,
        "only_source_ids": sorted(requested_ids),
        "controls": len(rows) - len(candidates),
        "note": "Crops are evidence for manual review; no automatic raster label recognition or accuracy claim.",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("toposheet", type=Path)
    parser.add_argument("--index", type=Path, default=Path("index_50k.geojson"))
    parser.add_argument("--sheet", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--context-m", type=float, default=600)
    parser.add_argument("--triangle-side", type=float, default=70)
    parser.add_argument(
        "--only-id",
        action="append",
        default=[],
        help="Restrict review to a source feature ID; repeatable",
    )
    args = parser.parse_args()
    audit(
        args.source,
        args.toposheet,
        args.index,
        args.sheet,
        args.output_dir,
        args.context_m,
        args.triangle_side,
        args.only_id,
    )
