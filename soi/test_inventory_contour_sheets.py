# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyarrow>=18", "pyproj>=3.6", "shapely>=2.1", "pytest>=8"]
# ///
import csv
import json

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import shapely
from shapely.geometry import LineString, MultiLineString, box, mapping

from inventory_contour_sheets import SheetInventory, inventory


def test_true_geometry_intersections_edges_touches_and_no_data():
    polygons = [box(i * 20, 0, i * 20 + 10, 10) for i in range(5)]
    collector = SheetInventory(polygons, [p.buffer(-1) for p in polygons])
    lines = [
        LineString([(1, 5), (9, 5)]),
        LineString([(20.1, 0.1), (29.9, 0.1)]),
        LineString([(0, 0), (0, 10)]),
        LineString([(35, 5), (40, 5)]),
        # Bbox covers the last sheet, but neither component intersects it.
        MultiLineString([[(75, 1), (76, 1)], [(95, 9), (96, 9)]]),
        LineString([(2, 2), (2, 2)]),
        None,
        LineString(),
    ]
    collector.add(lines)
    assert collector.intersects.tolist() == [2, 1, 1, 0, 0]
    assert collector.linework.tolist() == [2, 1, 0, 0, 0]
    assert collector.interior.tolist() == [1, 0, 0, 0, 0]
    features = [{"properties": {"id": str(i)}} for i in range(5)]
    assert [r["status"] for r in collector.records(features)] == [
        "interior",
        "edge_only",
        "touch_only",
        "no_data",
        "no_data",
    ]
    assert collector.statistics["zero_length_geometries"] == 1
    assert collector.statistics["null_geometries"] == 1
    assert collector.statistics["empty_geometries"] == 1


def test_one_source_row_can_cross_two_sheet_interiors():
    polygons = [box(0, 0, 10, 10), box(10, 0, 20, 10)]
    collector = SheetInventory(polygons, [p.buffer(-1) for p in polygons])
    collector.add([LineString([(5, 5), (15, 5)])])
    assert collector.interior.tolist() == [1, 1]
    assert collector.statistics["rows_intersecting_index"] == 1


def test_partitioned_streaming_inventory_and_metric_edge_band(tmp_path):
    polygons = [box(79 + i * 0.25, 20, 79 + (i + 1) * 0.25, 20.25) for i in range(4)]
    features = [
        {"type": "Feature", "properties": {"id": f"S{i}"}, "geometry": mapping(p)}
        for i, p in enumerate(polygons)
    ]
    index = tmp_path / "index.geojson"
    index.write_text(json.dumps({"type": "FeatureCollection", "features": features}))
    lines = [
        LineString([(79.1, 20.1), (79.2, 20.2)]),
        # ~1 m from the southern edge, wholly inside S1's two-metre band.
        LineString([(79.3, 20.000009), (79.4, 20.000009)]),
        # Only an endpoint reaches S2; most of the line is in S1.
        LineString([(79.45, 20.1), (79.5, 20.1)]),
    ]
    geo = {
        "version": "1.1.0",
        "primary_column": "geometry",
        "columns": {"geometry": {"encoding": "WKB", "geometry_types": ["LineString"]}},
    }
    extents = {}
    for n, pieces in enumerate([lines[:1], lines[1:]]):
        path = tmp_path / f"part{n}.parquet"
        table = pa.table(
            {"geometry": pa.array(shapely.to_wkb(pieces).tolist(), type=pa.binary())}
        )
        table = table.replace_schema_metadata({b"geo": json.dumps(geo).encode()})
        pq.write_table(table, path, row_group_size=1)
        extents[path.name] = {}
    metadata = tmp_path / "data.parquet.meta.json"
    metadata.write_text(json.dumps({"extents": extents}))
    output = tmp_path / "inventory"
    report = inventory([metadata], index, output, batch_size=1)
    assert report["source_statistics"]["source_rows"] == 3
    with (output / "sheet_coverage.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert [r["status"] for r in rows] == [
        "interior",
        "interior",
        "touch_only",
        "no_data",
    ]
    assert rows[1]["edge_only_source_rows"] == "1"
    assert rows[1]["interior_source_rows"] == "1"
    assert (output / "sheets_interior.txt").read_text() == "S0\nS1\n"
    assert (
        len(json.loads((output / "sheet_coverage.geojson").read_text())["features"])
        == 4
    )
    with pytest.raises(ValueError, match="empty output"):
        inventory([metadata], index, output)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
