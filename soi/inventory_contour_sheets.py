# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyarrow>=18", "pyproj>=3.6", "shapely>=2.1", "matplotlib>=3.9"]
# ///
"""Inventory SOI sheet coverage by streaming local GeoParquet geometries.

Counts refer to released source rows, including encoded height symbols. A row
can intersect several sheets. Presence does not establish complete DEM coverage.
"""

import argparse
import csv
import hashlib
import json
import math
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import shapely
from pyproj import CRS, Transformer
from shapely.geometry import shape


def parquet_paths(inputs):
    paths = []
    for value in inputs:
        path = Path(value).resolve()
        if path.name.endswith(".parquet.meta.json"):
            extents = json.loads(path.read_text())["extents"]
            if not extents:
                raise ValueError("Partition metadata has no files")
            for name in extents:
                if Path(name).name != name or not name.endswith(".parquet"):
                    raise ValueError(f"Invalid partition name: {name}")
                paths.append(path.parent / name)
        else:
            paths.append(path)
    if not paths or len(set(paths)) != len(paths):
        raise ValueError("Provide a nonempty set of distinct Parquet files")
    for path in paths:
        if not path.is_file() or path.suffix != ".parquet":
            raise ValueError(f"Missing Parquet: {path}")
    return paths


def geometry_column(parquet):
    metadata = json.loads((parquet.schema_arrow.metadata or {})[b"geo"])
    name = metadata["primary_column"]
    column = metadata["columns"][name]
    if column["encoding"] != "WKB" or column.get("edges", "planar") != "planar":
        raise ValueError("Expected planar WKB geometries")
    crs = column.get("crs", "OGC:CRS84")
    if crs is None or not CRS.from_user_input(crs).equals(
        CRS.from_epsg(4326), ignore_axis_order=True
    ):
        raise ValueError("Expected WGS84 longitude/latitude geometry")
    dtype = parquet.schema_arrow.field(name).type
    if not (pa.types.is_binary(dtype) or pa.types.is_large_binary(dtype)):
        raise ValueError("Expected a binary WKB geometry column")
    return name


def inset_sheets(polygons, margin_m):
    """Construct the inset in local WGS84 UTM, retaining curved sheet edges."""
    if not math.isfinite(margin_m) or margin_m < 0:
        raise ValueError("Edge margin must be finite and nonnegative")
    if margin_m == 0:
        return polygons.copy()
    transformers = {}
    result = []
    for polygon in polygons:
        centre = polygon.centroid
        zone = max(1, min(60, int((centre.x + 180) // 6) + 1))
        epsg = (32600 if centre.y >= 0 else 32700) + zone
        if epsg not in transformers:
            transformers[epsg] = (
                Transformer.from_crs(4326, epsg, always_xy=True),
                Transformer.from_crs(epsg, 4326, always_xy=True),
            )
        forward, backward = transformers[epsg]
        metric = shapely.transform(
            shapely.segmentize(polygon, 0.002), forward.transform, interleaved=False
        )
        inner = metric.buffer(-margin_m)
        if inner.is_empty:
            raise ValueError("Edge margin removes a sheet completely")
        result.append(
            shapely.transform(
                inner, backward.transform, interleaved=False
            ).intersection(polygon)
        )
    return np.asarray(result, dtype=object)


class SheetInventory:
    def __init__(self, polygons, interiors):
        self.polygons = np.asarray(polygons, dtype=object)
        self.interiors = np.asarray(interiors, dtype=object)
        shapely.prepare(self.interiors)
        self.tree = shapely.STRtree(self.polygons)
        self.intersects = np.zeros(len(polygons), dtype=np.int64)
        self.linework = np.zeros(len(polygons), dtype=np.int64)
        self.interior = np.zeros(len(polygons), dtype=np.int64)
        self.statistics = Counter()

    def add(self, geometries):
        geometries = np.asarray(geometries, dtype=object)
        self.statistics["source_rows"] += len(geometries)
        null = shapely.is_missing(geometries)
        empty = shapely.is_empty(geometries)
        self.statistics["null_geometries"] += int(null.sum())
        self.statistics["empty_geometries"] += int(empty.sum())
        types = shapely.get_type_id(geometries)
        nonline = ~null & ~empty & ~np.isin(types, [1, 5])
        if nonline.any():
            raise ValueError("Unexpected non-line geometry in contour input")
        lengths = shapely.length(geometries)
        unusable = ~null & ~empty & ~np.isfinite(lengths)
        if unusable.any():
            raise ValueError("Nonfinite contour geometry length")
        zero = ~null & ~empty & (lengths == 0)
        self.statistics["zero_length_geometries"] += int(zero.sum())
        geometries = geometries[~null & ~empty & ~zero]
        self.statistics["usable_line_rows"] += len(geometries)
        pairs = self.tree.query(geometries, predicate="intersects")
        source, sheet = pairs
        self.statistics["rows_intersecting_index"] += len(np.unique(source))
        self.statistics["rows_outside_index"] += len(geometries) - len(
            np.unique(source)
        )
        if not len(source):
            return
        np.add.at(self.intersects, sheet, 1)
        lines = geometries[source]
        # Most rows lie wholly inside one sheet. Prepared containment avoids
        # constructing intersection matrices against every densified edge.
        inside = shapely.contains_properly(self.interiors[sheet], lines)
        np.add.at(self.linework, sheet[inside], 1)
        np.add.at(self.interior, sheet[inside], 1)
        lines, sheet = lines[~inside], sheet[~inside]
        if not len(sheet):
            return
        # DE-9IM dimensions distinguish a positive-length piece from a point
        # touch, including linework coincident with the outer sheet boundary.
        relation = shapely.relate(lines, self.polygons[sheet])
        has_length = np.fromiter(
            (r[0] == "1" or r[1] == "1" for r in relation),
            dtype=bool,
            count=len(relation),
        )
        np.add.at(self.linework, sheet[has_length], 1)
        inside = shapely.relate_pattern(
            lines[has_length], self.interiors[sheet[has_length]], "1********"
        )
        np.add.at(self.interior, sheet[has_length][inside], 1)

    def records(self, features):
        records = []
        for i, feature in enumerate(features):
            any_rows, line_rows, inner_rows = map(
                int, (self.intersects[i], self.linework[i], self.interior[i])
            )
            status = (
                "interior"
                if inner_rows
                else "edge_only"
                if line_rows
                else "touch_only"
                if any_rows
                else "no_data"
            )
            records.append(
                {
                    "sheet": feature["properties"]["id"],
                    "osm_sheet": feature["properties"].get("OSM_SHEET_", ""),
                    "status": status,
                    "has_interior_data": bool(inner_rows),
                    "intersecting_source_rows": any_rows,
                    "linework_source_rows": line_rows,
                    "interior_source_rows": inner_rows,
                    "edge_only_source_rows": line_rows - inner_rows,
                    "point_touch_source_rows": any_rows - line_rows,
                }
            )
        return records


def plot_inventory(features, records, output, margin_m):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Patch, Polygon

    colours = {
        "interior": "#238b45",
        "edge_only": "#f0a13a",
        "touch_only": "#b484c4",
        "no_data": "#e2e5e9",
    }
    counts = Counter(r["status"] for r in records)
    fig, ax = plt.subplots(figsize=(12, 12), layout="constrained")
    fig.get_layout_engine().set(rect=(0, 0.04, 1, 0.96))
    patches = [Polygon(np.asarray(f["geometry"]["coordinates"][0])) for f in features]
    ax.add_collection(
        PatchCollection(
            patches,
            facecolors=[colours[r["status"]] for r in records],
            edgecolors="#69717b",
            linewidths=0.13,
        )
    )
    ax.autoscale()
    ax.set_aspect(1 / math.cos(math.radians(23)))
    ax.set(
        xlabel="Longitude (WGS84)",
        ylabel="Latitude (WGS84)",
        title="Released SOI contour data by 1:50,000 sheet",
    )
    labels = {
        "interior": f"Interior linework beyond {margin_m:g} m edge band",
        "edge_only": "Edge-band linework only",
        "touch_only": "Point contact only",
        "no_data": "No intersecting linework",
    }
    ax.legend(
        handles=[
            Patch(facecolor=colours[s], label=f"{labels[s]}: {counts[s]:,}")
            for s in colours
        ],
        loc="lower left",
        fontsize=8,
    )
    fig.text(
        0.5,
        0.015,
        "Presence of raw contour features, including height symbols; not a completeness or DEM-quality assessment.",
        ha="center",
        fontsize=8,
    )
    fig.savefig(output, dpi=180)
    plt.close(fig)


def inventory(
    inputs,
    index,
    output_dir,
    *,
    margin_m=2.0,
    batch_size=8192,
    plot=False,
    progress=print,
):
    if not isinstance(batch_size, int) or batch_size <= 0:
        raise ValueError("Batch size must be a positive integer")
    output_dir = Path(output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("Choose a new, empty output directory")
    paths = parquet_paths(inputs)
    features = json.loads(Path(index).read_text())["features"]
    ids = [f["properties"]["id"] for f in features]
    if not features or len(set(ids)) != len(ids):
        raise ValueError("Sheet IDs must be nonempty and unique")
    polygons = np.asarray([shape(f["geometry"]) for f in features], dtype=object)
    if any(g.geom_type != "Polygon" or g.is_empty or not g.is_valid for g in polygons):
        raise ValueError("Expected valid sheet polygons")
    started = time.monotonic()
    collector = SheetInventory(polygons, inset_sheets(polygons, margin_m))
    sources = []
    for path in paths:
        before = path.stat()
        with pq.ParquetFile(path) as parquet:
            column = geometry_column(parquet)
            total = parquet.metadata.num_rows
            progress(f"Reading {path.name}: {total:,} rows", flush=True)
            seen, last_update = 0, time.monotonic()
            for batch in parquet.iter_batches(
                batch_size=batch_size, columns=[column], use_threads=False
            ):
                collector.add(
                    shapely.from_wkb(batch.column(0).to_pylist(), on_invalid="raise")
                )
                seen += batch.num_rows
                if time.monotonic() - last_update >= 15:
                    progress(
                        f"{path.name}: {seen:,}/{total:,} rows; {np.count_nonzero(collector.interior):,} sheets with interior linework so far",
                        flush=True,
                    )
                    last_update = time.monotonic()
            if seen != total:
                raise ValueError(f"Unexpected row count in {path}")
        after = path.stat()
        if (before.st_size, before.st_mtime_ns, before.st_ino) != (
            after.st_size,
            after.st_mtime_ns,
            after.st_ino,
        ):
            raise ValueError(f"Input changed during scan: {path}")
        sources.append(
            {
                "path": str(path),
                "bytes": after.st_size,
                "mtime_ns": after.st_mtime_ns,
                "rows": seen,
            }
        )
    records = collector.records(features)
    output_dir.mkdir(parents=True, exist_ok=True)
    ordered = sorted(records, key=lambda r: r["sheet"])
    with (output_dir / "sheet_coverage.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(ordered)
    mapped = [
        dict(f, properties=dict(f["properties"], **r))
        for f, r in zip(features, records)
    ]
    (output_dir / "sheet_coverage.geojson").write_text(
        json.dumps({"type": "FeatureCollection", "features": mapped}) + "\n"
    )
    for status in ("interior", "edge_only", "touch_only", "no_data"):
        (output_dir / f"sheets_{status}.txt").write_text(
            "".join(r["sheet"] + "\n" for r in ordered if r["status"] == status)
        )
    report = {
        "index": str(index),
        "index_sha256": hashlib.sha256(Path(index).read_bytes()).hexdigest(),
        "sheet_grid_crs": "EPSG:4326",
        "edge_margin_m": margin_m,
        "sheet_count": len(features),
        "status_counts": dict(Counter(r["status"] for r in records)),
        "source_statistics": dict(collector.statistics),
        "sources": sources,
        "elapsed_seconds": time.monotonic() - started,
        "batch_size": batch_size,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "definitions": {
            "interior": "Positive-length raw linework intersects the sheet interior beyond the configured edge band.",
            "edge_only": "Positive-length raw linework exists, but only in the configured edge band or on the boundary.",
            "touch_only": "Source linework touches the sheet at points only; no positive-length intersection.",
            "no_data": "No usable source line geometry intersects the sheet.",
        },
        "limitations": [
            "Counts refer to released source rows, not deduplicated physical contours; one row can count in several sheets.",
            "Small circular/triangular height symbols remain included as raw linework.",
            "Presence does not establish complete coverage, usable height attributes, consistent elevations or DEM quality.",
        ],
    }
    (output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    if plot:
        plot_inventory(features, records, output_dir / "sheet_coverage.png", margin_m)
    progress(
        json.dumps(
            {
                k: report[k]
                for k in (
                    "sheet_count",
                    "status_counts",
                    "source_statistics",
                    "elapsed_seconds",
                )
            },
            indent=2,
        ),
        flush=True,
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument(
        "--index", type=Path, default=Path(__file__).with_name("index_50k.geojson")
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--edge-margin-m", type=float, default=2)
    parser.add_argument("--batch-size", type=int, default=8192)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    try:
        inventory(
            args.inputs,
            args.index,
            args.output_dir,
            margin_m=args.edge_margin_m,
            batch_size=args.batch_size,
            plot=args.plot,
        )
    except (ValueError, KeyError, OSError, shapely.errors.GEOSException) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
