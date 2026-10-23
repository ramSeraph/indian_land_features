# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "shapely>=2.1", "matplotlib>=3.9"]
# ///
"""Trace DEM conflict locations back to unchanged source contour linework.

Reads cached extracts only. Produces a JSON audit, source GeoJSONL, and plots;
never changes heights, clips contours, or resolves ambiguous constraints.
"""

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from shapely import STRtree, normalize, segmentize
from shapely import transform as transform_coordinates
from shapely.geometry import MultiLineString, Point, box, mapping, shape

TOLERANCE_M = 0.002  # Covers the DEM's millimetre coordinate rounding.


def read_features(path):
    with Path(path).open() as stream:
        for line in stream:
            if line.strip():
                yield json.loads(line.strip().lstrip("\x1e"))


def transform(function, geometry):
    return transform_coordinates(geometry, function, interleaved=False)


def load_lines(path, project):
    records = []
    for feature in read_features(path):
        geometry = shape(feature["geometry"])
        components = (
            geometry.geoms if geometry.geom_type == "MultiLineString" else [geometry]
        )
        for part_index, line in enumerate(components):
            if line.geom_type != "LineString":
                continue
            projected = transform(project, line)
            props = feature.get("properties") or {}
            records.append(
                {
                    "source_id": props.get("id", feature.get("id")),
                    "part_index": part_index,
                    "elevation_m": float(props["VALUE"]),
                    "source_properties": props,
                    "line": line,
                    "projected": projected,
                }
            )
    return records


def describe_line(record):
    line, projected = record["line"], record["projected"]
    return {
        key: record[key]
        for key in ("source_id", "part_index", "elevation_m", "source_properties")
    } | {
        "length_m": projected.length,
        "vertex_count": len(line.coords),
        "is_ring": line.is_ring,
        "is_simple": line.is_simple,
        "bounds": list(line.bounds),
        "endpoints": [list(line.coords[0]), list(line.coords[-1])],
        "geometry_sha256": hashlib.sha256(normalize(line).wkb).hexdigest(),
    }


def contact_kind(line, point):
    vertices = np.asarray(line.coords)[:, :2]
    distance = np.linalg.norm(vertices - [point.x, point.y], axis=1)
    if min(distance[0], distance[-1]) <= TOLERANCE_M:
        return "endpoint"
    if distance.min() <= TOLERANCE_M:
        return "interior_vertex"
    return "segment_interior"


def endpoint_neighbors(index, records, tree):
    record = records[index]
    result = []
    for xy in [record["projected"].coords[0], record["projected"].coords[-1]]:
        point = Point(xy)
        neighbors = []
        for other_index in tree.query(point.buffer(TOLERANCE_M)):
            other_index = int(other_index)
            other = records[other_index]
            if other_index == index or other["projected"].length == 0:
                continue
            if other["projected"].distance(point) <= TOLERANCE_M:
                neighbors.append(
                    {
                        "source_id": other["source_id"],
                        "part_index": other["part_index"],
                        "elevation_m": other["elevation_m"],
                        "length_m": other["projected"].length,
                        "contact": contact_kind(other["projected"], point),
                    }
                )
        result.append(neighbors)
    return result


def candidate_grids(index_path, sheet, project):
    """Test cutting-grid interpretations; leave the contour coordinates alone."""
    index = json.loads(index_path.read_text())
    matches = [f for f in index["features"] if f["properties"].get("id") == sheet]
    if len(matches) != 1:
        raise ValueError(f"Expected one sheet polygon for {sheet}")
    nominal = shape(matches[0]["geometry"])
    left, bottom, right, top = nominal.bounds
    # The 7.5-minute lines divide a standard 15-minute 50k sheet into quarters.
    internal = MultiLineString(
        [
            [(left, (bottom + top) / 2), (right, (bottom + top) / 2)],
            [((left + right) / 2, bottom), ((left + right) / 2, top)],
        ]
    )
    candidates = []
    for epsg in (4326, 4145, 4146):
        conversion = Transformer.from_crs(
            epsg, 4326, always_xy=True, allow_ballpark=False
        )
        polygon = transform(conversion.transform, segmentize(nominal, 0.001))
        midlines = transform(conversion.transform, segmentize(internal, 0.001))
        outer_projected = transform(project, polygon)
        internal_projected = transform(project, midlines)
        operation = conversion.get_last_used_operation() if epsg != 4326 else conversion
        candidates.append(
            {
                "source_epsg": epsg,
                "operation": operation.description,
                "operation_accuracy_m": operation.accuracy,
                "operation_definition": operation.definition,
                "polygon": polygon,
                "internal": midlines,
                "projected_polygon": outer_projected,
                "projected_internal": internal_projected,
            }
        )
    return candidates


def audit(raw_path, classified_path, conflicts_path, output, sheet, index_path):
    conflicts = list(read_features(conflicts_path))
    if not conflicts:
        raise ValueError("No conflict locations supplied")
    centre = np.mean([f["geometry"]["coordinates"][:2] for f in conflicts], axis=0)
    epsg = (32600 if centre[1] >= 0 else 32700) + int((centre[0] + 180) // 6) + 1
    project = Transformer.from_crs(4326, epsg, always_xy=True).transform
    grids = candidate_grids(index_path, sheet, project)
    records = load_lines(raw_path, project)
    active = load_lines(classified_path, project)
    tree = STRtree([r["projected"] for r in records])
    active_tree = STRtree([r["projected"] for r in active])
    locations, pairs, involved = [], {}, set()
    for number, feature in enumerate(conflicts, 1):
        point = transform(project, shape(feature["geometry"]))
        active_hits = [
            int(i)
            for i in active_tree.query(point.buffer(TOLERANCE_M))
            if active[int(i)]["projected"].distance(point) <= TOLERANCE_M
        ]
        active_keys = {
            (active[i]["source_id"], active[i]["elevation_m"]) for i in active_hits
        }
        hits = [
            int(i)
            for i in tree.query(point.buffer(TOLERANCE_M))
            if records[int(i)]["projected"].length > 0
            and records[int(i)]["projected"].distance(point) <= TOLERANCE_M
            and (records[int(i)]["source_id"], records[int(i)]["elevation_m"])
            in active_keys
        ]
        hits.sort()
        row = {
            "conflict_id": f"C{number:03}",
            "coordinates": feature["geometry"]["coordinates"],
            **feature["properties"],
            "source_contacts": [
                describe_line(records[i])
                | {
                    "contact": contact_kind(records[i]["projected"], point),
                    "distance_m": records[i]["projected"].distance(point),
                }
                for i in hits
            ],
            "pairs": [],
            "candidate_grids": [
                {
                    "source_epsg": g["source_epsg"],
                    "outer_distance_m": g["projected_polygon"].boundary.distance(point),
                    "quarter_distance_m": g["projected_internal"].distance(point),
                    "inside_sheet": g["projected_polygon"].covers(point),
                }
                for g in grids
            ],
        }
        kinds = [c["contact"] for c in row["source_contacts"]]
        row["junction_type"] = (
            "two_line_height_change"
            if kinds == ["endpoint", "endpoint"]
            else "four_endpoint_junction"
            if kinds == ["endpoint"] * 4
            else "endpoint_meets_continuing_line"
            if sorted(kinds) == ["endpoint", "interior_vertex"]
            else "other"
        )
        involved.update(hits)
        for a, b in itertools.combinations(hits, 2):
            first, second = records[a], records[b]
            if abs(first["elevation_m"] - second["elevation_m"]) <= 0.01:
                continue
            key = (a, b)
            if key not in pairs:
                # Intersect the original coordinates, before any projection or clipping.
                intersection = first["line"].intersection(second["line"])
                overlap = transform(project, intersection).length
                equal = first["line"].equals(second["line"])
                covers = first["line"].covers(second["line"]) or second["line"].covers(
                    first["line"]
                )
                kinds = sorted(
                    contact_kind(r["projected"], point) for r in (first, second)
                )
                category = (
                    "identical_linework_different_heights"
                    if equal
                    else "overlapping_linework_different_heights"
                    if overlap > 0.001
                    else "_to_".join(kinds)
                )
                pairs[key] = {
                    "pair_id": f"P{len(pairs) + 1:03}",
                    "category": category,
                    "lines": [describe_line(first), describe_line(second)],
                    "source_intersection_type": intersection.geom_type,
                    "source_intersection_empty": intersection.is_empty,
                    "overlap_length_m": overlap,
                    "one_line_covers_other": covers,
                    "elevation_difference_m": abs(
                        first["elevation_m"] - second["elevation_m"]
                    ),
                    "conflicts": [],
                    "record_indices": [a, b],
                }
            pairs[key]["conflicts"].append(row["conflict_id"])
            row["pairs"].append(pairs[key]["pair_id"])
        locations.append(row)
    output.mkdir(parents=True, exist_ok=True)
    source_records = []
    for i in sorted(involved):
        record = records[i]
        props = describe_line(record) | {
            "endpoint_neighbors": endpoint_neighbors(i, records, tree)
        }
        source_records.append(
            {
                "type": "Feature",
                "properties": props,
                "geometry": mapping(record["line"]),
            }
        )
    with (output / "involved_source_lines.geojsonl").open("w") as stream:
        for row in source_records:
            stream.write(json.dumps(row) + "\n")
    report = {
        "sheet": sheet,
        "projected_epsg": epsg,
        "contact_tolerance_m": TOLERANCE_M,
        "inputs": [
            {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in (raw_path, classified_path, conflicts_path)
        ],
        "summary": {
            "conflict_locations": len(locations),
            "unmatched_locations": [
                r["conflict_id"] for r in locations if not r["pairs"]
            ],
            "source_line_pairs": len(pairs),
            "involved_source_line_parts": len(involved),
            "pair_categories": dict(Counter(p["category"] for p in pairs.values())),
            "location_junction_types": dict(
                Counter(r["junction_type"] for r in locations)
            ),
            "all_pairs_intersect_in_original_source_coordinates": all(
                not p["source_intersection_empty"] for p in pairs.values()
            ),
            "candidate_grid_counts": [
                {
                    "source_epsg": g["source_epsg"],
                    "outer_within_2m": sum(
                        r["candidate_grids"][i]["outer_distance_m"] <= 2
                        for r in locations
                    ),
                    "outer_within_25m": sum(
                        r["candidate_grids"][i]["outer_distance_m"] <= 25
                        for r in locations
                    ),
                    "outer_or_quarter_within_25m": sum(
                        min(
                            r["candidate_grids"][i]["outer_distance_m"],
                            r["candidate_grids"][i]["quarter_distance_m"],
                        )
                        <= 25
                        for r in locations
                    ),
                    "existing_locations_inside_candidate_sheet": sum(
                        r["candidate_grids"][i]["inside_sheet"] for r in locations
                    ),
                }
                for i, g in enumerate(grids)
            ],
        },
        "candidate_grid_operations": [
            {
                k: g[k]
                for k in (
                    "source_epsg",
                    "operation",
                    "operation_accuracy_m",
                    "operation_definition",
                )
            }
            for g in grids
        ],
        "pairs": list(pairs.values()),
        "locations": locations,
        "limitations": [
            "Audits the supplied DEM sample conflicts, not all possible unsampled contour crossings.",
            "Raw source matching excludes zero-length parts and requires a contacting retained contour of the same ID and height.",
            "Shape alone does not establish which conflicting elevation is correct or the source's original cutting datum.",
            "Candidate grids reinterpret index coordinates only; their original datum is unverified. Inside counts concern existing conflict locations, not a complete re-extraction or repaired DEM.",
        ],
    }
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    grid_features = []
    for g in grids:
        for kind in ("polygon", "internal"):
            grid_features.append(
                {
                    "type": "Feature",
                    "properties": {"source_epsg": g["source_epsg"], "grid_part": kind},
                    "geometry": mapping(g[kind]),
                }
            )
    (output / "candidate_grids.geojson").write_text(
        json.dumps({"type": "FeatureCollection", "features": grid_features}) + "\n"
    )
    plot_pairs(report, records, tree, output, project)
    plot_overview(report, records, grids, output)
    print(
        json.dumps(
            {"sheet": sheet, **report["summary"], "output": str(output)}, indent=2
        ),
        flush=True,
    )


def plot_pairs(report, records, tree, output, project):
    locations = {r["conflict_id"]: r for r in report["locations"]}
    palette = ["#d14900", "#0066cc"]
    for page in range(0, len(report["pairs"]), 6):
        fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)
        for ax, pair in zip(axes.flat, report["pairs"][page : page + 6]):
            points = np.array(
                [project(*locations[c]["coordinates"][:2]) for c in pair["conflicts"]]
            )
            centre = points.mean(axis=0)
            radius = max(50, min(500, np.max(np.abs(points - centre)) + 30))
            window = box(
                centre[0] - radius,
                centre[1] - radius,
                centre[0] + radius,
                centre[1] + radius,
            )
            for i in tree.query(window):
                line = records[int(i)]["projected"]
                if line.length == 0:
                    continue
                xy = np.asarray(line.coords)[:, :2] - centre
                ax.plot(xy[:, 0], xy[:, 1], color="0.78", lw=0.7, zorder=1)
            for color, index in zip(palette, pair["record_indices"]):
                record = records[index]
                xy = np.asarray(record["projected"].coords)[:, :2] - centre
                label = f"{record['source_id'].split('.')[-1]} / {record['elevation_m']:g} m"
                ax.plot(
                    xy[:, 0],
                    xy[:, 1],
                    color=color,
                    lw=1.5,
                    marker=".",
                    ms=2,
                    label=label,
                    zorder=3,
                )
                ax.scatter(
                    xy[[0, -1], 0],
                    xy[[0, -1], 1],
                    facecolors="none",
                    edgecolors=color,
                    s=45,
                    zorder=4,
                )
            local = points - centre
            ax.scatter(
                local[:, 0],
                local[:, 1],
                marker="x",
                color="black",
                s=45,
                linewidths=1.2,
                zorder=5,
            )
            ax.set(
                xlim=(-radius, radius),
                ylim=(-radius, radius),
                aspect="equal",
                xlabel="East from conflict centre (m)",
                ylabel="North (m)",
            )
            lengths = " / ".join(f"{r['length_m']:.2f}" for r in pair["lines"])
            title = f"{pair['pair_id']} · {', '.join(pair['conflicts'])}\n{pair['category'].replace('_', ' ')}\nLengths {lengths} m"
            ax.set_title(title, fontsize=8)
            ax.legend(fontsize=7, loc="best")
            ax.grid(alpha=0.2)
        for ax in axes.flat[len(report["pairs"][page : page + 6]) :]:
            ax.set_visible(False)
        fig.suptitle(
            f"{report['sheet']}: original contour lines — circles mark endpoints, × marks conflicts",
            fontsize=12,
        )
        fig.savefig(output / f"line-pairs-{page // 6 + 1:02}.png", dpi=150)
        plt.close(fig)


def plot_overview(report, records, grids, output):
    fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)
    for record in records:
        xy = np.asarray(record["line"].coords)
        ax.plot(xy[:, 0], xy[:, 1], color="0.8", lw=0.25)
    for grid, color, name in zip(
        grids[:2],
        ("#0066cc", "#b3261e"),
        ("WGS84 index", "Index interpreted as Kalianpur 1962"),
    ):
        for component in [grid["polygon"].exterior, *grid["internal"].geoms]:
            xy = np.asarray(component.coords)
            ax.plot(xy[:, 0], xy[:, 1], color=color, lw=0.8, ls="--")
        ax.plot([], [], color=color, lw=1, ls="--", label=name)
    xy = np.array([r["coordinates"][:2] for r in report["locations"]])
    ax.scatter(
        xy[:, 0],
        xy[:, 1],
        color="black",
        marker="x",
        s=25,
        label="Conflicting heights",
        zorder=5,
    )
    for row in report["locations"]:
        ax.annotate(
            row["conflict_id"][1:],
            row["coordinates"][:2],
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=6,
        )
    left, bottom, right, top = grids[0]["polygon"].bounds
    ax.set(
        xlim=(left - 0.005, right + 0.005),
        ylim=(bottom - 0.005, top + 0.005),
        aspect=1 / np.cos(np.deg2rad((bottom + top) / 2)),
        xlabel="Longitude (WGS84)",
        ylabel="Latitude (WGS84)",
        title=f"{report['sheet']}: original contours and conflict locations\nCandidate cutting grids; contour coordinates unchanged",
    )
    ax.legend(loc="best", fontsize=8)
    ax.ticklabel_format(useOffset=False)
    fig.savefig(output / "overview.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("classified", type=Path)
    parser.add_argument("conflicts", type=Path)
    parser.add_argument("--sheet", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--index", type=Path, default=Path(__file__).with_name("index_50k.geojson")
    )
    args = parser.parse_args()
    audit(
        args.raw,
        args.classified,
        args.conflicts,
        args.output_dir,
        args.sheet,
        args.index,
    )


if __name__ == "__main__":
    main()
