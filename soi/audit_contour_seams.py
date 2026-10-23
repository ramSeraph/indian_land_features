# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "rasterio>=1.4", "scipy>=1.14", "shapely>=2.1", "matplotlib>=3.9"]
# ///
"""Compare unchanged source heights at exact joins near a sheet's outer edges.

This audits the cached extract, including complete lines from neighbouring
sheets. It does not correct heights or assume that clipping repairs a seam.
"""

import argparse
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer
from shapely import segmentize
from shapely.geometry import LineString, Point, mapping, shape

from contours_to_dem import geometry_hash, get_sheet, parts, read_features, transform


def audit_seams(inputs, index, sheet_id, tolerance_m=2.0):
    sheet = get_sheet(index, sheet_id)
    zone = min(60, int((sheet.centroid.x + 180) // 6) + 1)
    epsg = (32600 if sheet.centroid.y >= 0 else 32700) + zone
    project = Transformer.from_crs(4326, epsg, always_xy=True).transform
    projected_sheet = transform(project, segmentize(sheet, 0.002))
    interior = projected_sheet.buffer(-tolerance_m)
    w, s, e, n = sheet.bounds
    edges = {
        "west": LineString([(w, s), (w, n)]),
        "east": LineString([(e, s), (e, n)]),
        "south": LineString([(w, s), (e, s)]),
        "north": LineString([(w, n), (e, n)]),
    }
    projected_edges = {
        name: transform(project, segmentize(line, 0.002))
        for name, line in edges.items()
    }
    # Neighbour names are from the nominal index, not an assertion of source
    # ownership: the released rows have no source sheet identifier.
    index_features = json.loads(Path(index).read_text())["features"]
    neighbours = {}
    probes = {
        "west": Point(w - 0.001, (s + n) / 2),
        "east": Point(e + 0.001, (s + n) / 2),
        "south": Point((w + e) / 2, s - 0.001),
        "north": Point((w + e) / 2, n + 0.001),
    }
    for name, probe in probes.items():
        neighbours[name] = [
            f["properties"].get("id")
            for f in index_features
            if shape(f["geometry"]).covers(probe)
        ]
    endpoints, records, seen = defaultdict(list), [], set()
    for path in inputs:
        for feature in read_features(path):
            props = feature.get("properties") or {}
            elevation = float(props["VALUE"])
            for part_index, line in enumerate(
                parts(shape(feature["geometry"]), "LineString")
            ):
                digest = geometry_hash(line)
                key = digest, elevation
                if key in seen or line.length == 0 or line.is_closed:
                    continue
                seen.add(key)
                projected = transform(project, segmentize(line, 0.002))
                locations = [tuple(line.coords[0][:2]), tuple(line.coords[-1][:2])]
                nearby = [
                    xy
                    for xy in locations
                    if projected_sheet.boundary.distance(Point(project(*xy)))
                    <= tolerance_m
                ]
                if not nearby:
                    continue
                clipped = line.intersection(sheet)
                retained = transform(project, segmentize(clipped, 0.002))
                inside_fraction = retained.length / projected.length
                record = {
                    "source_id": props.get("id", feature.get("id")),
                    "source_part_index": part_index,
                    "source_geometry_sha256": digest,
                    "elevation_m": elevation,
                    "source_length_m": projected.length,
                    "retained_length_m": retained.length,
                    "retained_fraction": inside_fraction,
                    "majority_side": "inside" if inside_fraction >= 0.5 else "outside",
                    "outside_edge_sliver": bool(
                        0 < inside_fraction < 0.5 and not retained.intersects(interior)
                    ),
                    "source_properties": props,
                    "geometry": mapping(line),
                }
                record_index = len(records)
                records.append(record)
                for xy in nearby:
                    endpoints[xy].append(record_index)
    joins, unmatched, ambiguous = [], [], []
    for xy, indices in sorted(endpoints.items()):
        point = Point(project(*xy))
        edge = min(projected_edges, key=lambda k: projected_edges[k].distance(point))
        inside = [i for i in indices if records[i]["majority_side"] == "inside"]
        outside = [i for i in indices if records[i]["majority_side"] == "outside"]
        base = {
            "coordinates": list(xy),
            "edge": edge,
            "nominal_neighbours": neighbours[edge],
            "distance_to_edge_m": projected_edges[edge].distance(point),
        }
        if not inside or not outside:
            unmatched.append(
                base | {"source_ids": [records[i]["source_id"] for i in indices]}
            )
            continue
        if len(inside) != 1 or len(outside) != 1:
            ambiguous.append(
                base | {"source_ids": [records[i]["source_id"] for i in indices]}
            )
            continue
        for i, j in itertools.product(inside, outside):
            delta = records[i]["elevation_m"] - records[j]["elevation_m"]
            joins.append(
                base
                | {
                    "join_id": f"J{len(joins) + 1:04}",
                    "inside": {k: v for k, v in records[i].items() if k != "geometry"},
                    "outside": {k: v for k, v in records[j].items() if k != "geometry"},
                    "inside_minus_outside_m": delta,
                    "height_conflict": abs(delta) > 0.01,
                    "record_indices": [i, j],
                }
            )
    summary = {}
    for edge in edges:
        rows = [row for row in joins if row["edge"] == edge]
        summary[edge] = {
            "nominal_neighbours": neighbours[edge],
            "matched_locations": len(rows),
            "conflicting_locations": sum(row["height_conflict"] for row in rows),
            "inside_minus_outside_counts": dict(
                sorted(Counter(row["inside_minus_outside_m"] for row in rows).items())
            ),
        }
    return {
        "sheet": sheet_id,
        "projected_epsg": epsg,
        "edge_tolerance_m": tolerance_m,
        "matched_locations": len(joins),
        "conflicting_locations": sum(row["height_conflict"] for row in joins),
        "inputs": [
            {"path": str(p), "sha256": hashlib.sha256(Path(p).read_bytes()).hexdigest()}
            for p in [*inputs, index]
        ],
        "summary": summary,
        "joins": joins,
        "unmatched_endpoint_groups": unmatched,
        "ambiguous_endpoint_groups": ambiguous,
        "limitations": [
            "Only exact source endpoint joins within the edge tolerance are compared.",
            "The sheet extract is not a complete extract of its neighbours; unmatched or absent joins are not evidence of agreement.",
            "Inside/outside attribution uses the majority of each uncut line's length, not a source sheet identifier.",
            "Excluding an outside sliver from one sheet does not reconcile the two source elevations for a mosaic.",
            "No heights or source coordinates are changed by this audit.",
        ],
    }, records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--sheet", required=True)
    parser.add_argument(
        "--index", type=Path, default=Path(__file__).with_name("index_50k.geojson")
    )
    parser.add_argument("--edge-tolerance", type=float, default=2.0)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not 0 < args.edge_tolerance <= 25:
        parser.error("Edge tolerance must be positive and no greater than 25 m")
    report, records = audit_seams(
        args.inputs, args.index, args.sheet, args.edge_tolerance
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    with (args.output_dir / "joins.geojsonl").open("w") as stream:
        for row in report["joins"]:
            feature = {
                "type": "Feature",
                "geometry": mapping(Point(row["coordinates"])),
                "properties": {k: v for k, v in row.items() if k != "record_indices"},
            }
            stream.write(json.dumps(feature) + "\n")
    with (args.output_dir / "source_lines.geojsonl").open("w") as stream:
        for row in records:
            stream.write(
                json.dumps(
                    {
                        "type": "Feature",
                        "geometry": row["geometry"],
                        "properties": {k: v for k, v in row.items() if k != "geometry"},
                    }
                )
                + "\n"
            )
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    differences = [row["inside_minus_outside_m"] for row in report["joins"]] or [0]
    limits = min(0, min(differences)) - 5, max(0, max(differences)) + 5
    for ax, edge in zip(axes.flat, ["west", "east", "south", "north"]):
        rows = [row for row in report["joins"] if row["edge"] == edge]
        coordinate = 1 if edge in {"west", "east"} else 0
        ax.axhline(0, color="#999999", lw=0.8)
        for row in rows:
            x, y = row["coordinates"][coordinate], row["inside_minus_outside_m"]
            ax.plot(
                x,
                y,
                "o",
                color="#c02020" if row["height_conflict"] else "#008050",
                ms=4,
            )
        ax.set_title(
            f"{edge}: {len(rows)} exact joins, {sum(r['height_conflict'] for r in rows)} mismatches"
        )
        ax.set_xlabel("Latitude (degrees)" if coordinate else "Longitude (degrees)")
        ax.set_ylabel("Inside minus outside height (m)")
        ax.set_ylim(*limits)
        ax.ticklabel_format(useOffset=False, style="plain")
        if not rows:
            ax.text(
                0.5,
                0.5,
                "No paired joins in cached extract",
                ha="center",
                transform=ax.transAxes,
            )
    fig.suptitle(f"{args.sheet}: source contour heights across nominal sheet edges")
    fig.savefig(args.output_dir / "height-jumps.png", dpi=150)
    plt.close(fig)
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
