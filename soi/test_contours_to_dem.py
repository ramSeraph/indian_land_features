# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "rasterio>=1.4", "scipy>=1.14", "shapely>=2.1", "pytest>=8"]
# ///
import json
from functools import partial

import numpy as np
import pytest
import rasterio
from pyproj import Transformer
from rasterio.transform import from_origin
from scipy.interpolate import LinearNDInterpolator
from shapely import segmentize
from shapely.affinity import rotate, scale
from shapely.geometry import LineString, MultiLineString, Point, Polygon, box, mapping

from contours_to_dem import (
    build_dem as build_sheet_dem,
    flatten_water,
    geometry_hash,
    prepare_constraints,
    roundness,
    spot_ring,
    spot_triangle,
    transform,
    unique_samples,
)

# Preserve the historical TIN regression suite; the default Imagico engine has
# its own integration checks in test_imagico_dem.py.
build_dem = partial(build_sheet_dem, interpolation="linear")


def feature(geometry, elevation):
    return {
        "type": "Feature",
        "properties": {"VALUE": elevation},
        "geometry": mapping(geometry),
    }


def test_edge_slivers_exclude_only_mostly_outside_shallow_parts(tmp_path):
    sheet = box(0, 0, 0.01, 0.01)
    project = lambda x, y: (x * 100000, y * 100000)
    source = tmp_path / "contours.geojsonl"
    sliver = LineString([(-0.001, 0.003), (0.00001, 0.003)])
    deep_crossing = LineString([(-0.001, 0.004), (0.00005, 0.004)])
    inside_edge = LineString([(0.00001, 0.006), (0.00001, 0.007)])
    # The circle is deliberately inside the 2 m edge band and must remain a point.
    circle = LineString(
        Point(0.00001, 0.008).buffer(0.000005, quad_segs=16).exterior.coords
    )
    source.write_text(
        "\n".join(
            json.dumps(feature(g, z))
            for g, z in [
                (sliver, 20),
                (deep_crossing, 40),
                (inside_edge, 60),
                (circle, 81),
            ]
        )
    )
    excluded = []
    contours, spots, _, counts, _ = prepare_constraints(
        [source],
        sheet,
        project,
        edge_sliver_m=2,
        excluded_contours=excluded,
    )
    assert [p["VALUE"] for _, p in contours] == [40, 60]
    assert [p["VALUE"] for _, p in spots] == [81]
    assert counts["excluded_outside_edge_sliver_parts"] == 1
    assert excluded[0][1]["source_geometry_sha256"] == geometry_hash(sliver)
    assert excluded[0][1]["exclusion_evidence"]["seam_reconciled"] is False
    original, _, _, _, _ = prepare_constraints([source], sheet, project)
    assert [p["VALUE"] for _, p in original] == [20, 40, 60]


def test_form_lines_require_uncut_t_junction_and_keep_other_half_interval_data(
    tmp_path,
):
    sheet = box(0, 0, 0.01, 0.01)
    project = lambda x, y: (x * 100000, y * 100000)
    regular = LineString([(0.004, 0.002), (0.004, 0.005), (0.004, 0.008)])
    branch = LineString([(0.004, 0.005), (0.002, 0.006)])
    endpoint_join = LineString([(0.004, 0.002), (0.003, 0.001)])
    isolated = LineString([(0.006, 0.004), (0.007, 0.004)])
    circle = LineString(
        Point(0.007, 0.007).buffer(0.00008, quad_segs=16).exterior.coords
    )
    source = tmp_path / "contours.geojsonl"
    source.write_text(
        "\n".join(
            json.dumps(feature(g, z))
            for g, z in [
                (regular, 100),
                (branch, 90),
                (endpoint_join, 90),
                (isolated, 90),
                (circle, 90),
            ]
        )
    )
    excluded = []
    contours, spots, _, counts, _ = prepare_constraints(
        [source],
        sheet,
        project,
        form_line_interval=20,
        excluded_contours=excluded,
    )
    assert len(contours) == 3
    assert len(spots) == 1 and spots[0][1]["VALUE"] == 90
    assert counts["excluded_form_line_candidate_parts"] == 1
    assert excluded[0][1]["source_geometry_sha256"] == geometry_hash(branch)
    assert excluded[0][1]["VALUE"] == 90
    # With a 10 m regular interval, 90 is an ordinary level, not a candidate.
    contours, _, _, _, _ = prepare_constraints(
        [source], sheet, project, form_line_interval=10
    )
    assert len(contours) == 4
    with pytest.raises(ValueError, match="interval"):
        prepare_constraints([source], sheet, project, form_line_interval=0)


@pytest.mark.parametrize("sparse", [False, True])
def test_automatic_form_line_audit_survives_dem_success_or_conflict(tmp_path, sparse):
    unproject = Transformer.from_crs(32644, 4326, always_xy=True).transform
    x, y = 540000, 2200000
    sheet = transform(unproject, box(x, y, x + 500, y + 500))
    index = tmp_path / "index.geojson"
    index.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"id": "test"},
                        "geometry": mapping(sheet),
                    }
                ],
            }
        )
    )
    positions = [(250, 100)] if sparse else [(50, 80), (250, 100), (450, 120)]
    lines = [
        (LineString([(x + dx, y - 2000), (x + dx, y + 250), (x + dx, y + 2500)]), z)
        for dx, z in positions
    ]
    lines.append((LineString([(x + 250, y + 250), (x + 150, y + 300)]), 90))
    if not sparse:
        lines.extend(
            [
                (LineString(Point(x + 350, y + 350).buffer(8).exterior.coords), 108),
                (
                    LineString(
                        [
                            (x + 70, y + 70),
                            (x + 140, y + 70),
                            (x + 105, y + 131),
                            (x + 70, y + 70),
                        ]
                    ),
                    89,
                ),
            ]
        )
    source = tmp_path / "source.geojsonl"
    source.write_text(
        "\n".join(json.dumps(feature(transform(unproject, g), z)) for g, z in lines)
    )
    output = tmp_path / "dem"
    if sparse:
        with pytest.raises(ValueError, match="conflicting elevations"):
            build_dem([source], index, "test", output)
        assert not (output / "dem.tif").exists()
    else:
        report = build_dem([source], index, "test", output)
        assert report["form_line_detection"]["mode"] == "auto"
        assert (output / "dem.tif").exists()
        points = [
            json.loads(s) for s in (output / "spots.geojsonl").read_text().splitlines()
        ]
        assert {
            (p["properties"]["symbol_type"], p["properties"]["VALUE"]) for p in points
        } == {("circle", 108), ("triangle", 89)}
    audit = [
        json.loads(s)
        for s in (output / "form_line_candidates.geojsonl").read_text().splitlines()
    ]
    assert len(audit) == 1
    assert audit[0]["properties"]["form_line_evidence"]["decision"] == (
        "retain_uncertain_interval" if sparse else "exclude_form_line_candidate"
    )


def test_ring_requires_closed_circle_and_accepts_both_symbol_sizes():
    for radius in (9, 11.5):
        ring = LineString(
            Point(500000, 2000000).buffer(radius, quad_segs=16).exterior.coords
        )
        point, fitted = spot_ring(ring)
        assert point.distance(Point(500000, 2000000)) < 1e-6
        assert fitted == pytest.approx(radius)
        assert spot_ring(LineString(list(ring.coords)[:-1])) is None
    assert spot_ring(LineString(Point(0, 0).buffer(30).exterior.coords)) is None


def test_dem_handles_cartographic_contact_and_keeps_height_point_inside_it(tmp_path):
    unproject = Transformer.from_crs(32644, 4326, always_xy=True).transform
    x, y = 540000, 2200000
    sheet = transform(unproject, box(x, y, x + 1000, y + 1000))
    index = tmp_path / "index.geojson"
    index.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"id": "test"},
                        "geometry": mapping(sheet),
                    }
                ],
            }
        )
    )
    source = tmp_path / "source.geojsonl"
    lines = [
        (LineString([(x, y), (x + 500, y + 500), (x + 1000, y)]), 80),
        (LineString([(x, y + 1000), (x + 500, y + 500), (x + 1000, y + 1000)]), 100),
        (LineString(Point(x + 500, y + 502).buffer(8).exterior.coords), 93),
    ]
    source.write_text(
        "\n".join(json.dumps(feature(transform(unproject, g), z)) for g, z in lines)
    )
    original = source.read_bytes()
    output = tmp_path / "dem"
    report = build_dem([source], index, "test", output)
    assert report["contour_contact_handling"]["contact_locations"] == 1
    assert report["contour_contact_handling"]["radius_m"] == 15
    assert report["counts"]["spot_heights"] == 1
    assert source.read_bytes() == original
    points = [
        json.loads(line)
        for line in (output / "spots.geojsonl").read_text().splitlines()
    ]
    assert points[0]["properties"]["VALUE"] == 93
    original_contours = [
        json.loads(line)
        for line in (output / "contours.geojsonl").read_text().splitlines()
    ]
    assert {f["properties"]["VALUE"] for f in original_contours} == {80, 100}
    assert all(len(f["geometry"]["coordinates"]) >= 3 for f in original_contours)
    with rasterio.open(output / "dem.tif") as raster:
        values = raster.read(1, masked=True)
        assert values.min() >= 80 and values.max() <= 100
    with pytest.raises(ValueError, match="conflicting elevations"):
        build_dem(
            [source], index, "test", tmp_path / "strict", contour_contact_radius_m=0
        )


def test_roundness_rejects_arc_closed_with_a_chord():
    angles = np.linspace(0, np.pi, 40)
    points = np.column_stack([9 * np.cos(angles), 9 * np.sin(angles)])
    closed_arc = LineString(np.vstack([points, points[0]]))
    # Every vertex lies on a 9 m circle, so radius variance alone misses this.
    assert closed_arc.is_ring
    assert roundness(closed_arc) < 0.8
    assert spot_ring(closed_arc) is None


def test_triangle_accepts_rotation_winding_and_extra_edge_samples():
    triangle = LineString([(100, 100), (165, 100), (132.5, 165), (100, 100)])
    for angle in (0, 37, 123):
        rotated = rotate(triangle, angle)
        for sampled in (rotated, segmentize(rotated, 5)):
            vertices = list(sampled.coords)[:-1]
            # Also start halfway along a side, rather than at a corner.
            vertices = vertices[1:] + vertices[:1]
            for ordered in (vertices, list(reversed(vertices))):
                line = LineString([*ordered, ordered[0]])
                centre, measurements = spot_triangle(line)
                assert centre.distance(Polygon(rotated).centroid) < 1e-8
                assert measurements["symbol_side_min_m"] == pytest.approx(65)
                assert measurements["symbol_edge_error_m"] < 1e-8
    noisy = LineString(
        [(0, 0), (32.5, -0.3), (65, 0), (49, 32.6), (32.5, 65), (16.1, 32.4), (0, 0)]
    )
    assert spot_triangle(noisy) is not None


def test_triangle_rejects_other_shapes_scales_and_nonstraight_sides():
    triangle = LineString([(0, 0), (65, 0), (32.5, 65), (0, 0)])
    bad_shapes = [
        LineString([(0, 0), (0, 0)]),
        LineString(list(triangle.coords)[:-1]),
        LineString([(0, 0), (65, 0), (32.5, 65), (0, 0.01)]),
        scale(triangle, xfact=0.4, yfact=0.4),
        scale(triangle, xfact=3, yfact=3),
        LineString([(0, 0), (70, 0), (35, 5), (0, 0)]),
        LineString([(0, 0), (65, 0), (65, 65), (0, 65), (0, 0)]),
        LineString([(0, 0), (65, 65), (0, 65), (65, 0), (0, 0)]),
        # Its convex hull is a correctly sized triangle, but one side is bent.
        LineString([(0, 0), (32.5, 5), (65, 0), (32.5, 65), (0, 0)]),
        LineString(Point(0, 0).buffer(30).exterior.coords),
    ]
    for line in bad_shapes:
        assert spot_triangle(line) is None
    assert spot_triangle(triangle, nominal_side_m=0) is None
    assert (
        spot_triangle(scale(triangle, xfact=3, yfact=3), nominal_side_m=210) is not None
    )


def test_real_56m5_triangle_is_detected_without_review_or_source_id(tmp_path):
    triangle = LineString(
        [
            (79.3352393, 19.8272999),
            (79.3355434, 19.8278893),
            (79.3358598, 19.8273057),
            (79.3352393, 19.8272999),
        ]
    )
    path = tmp_path / "source.geojsonl"
    path.write_text(json.dumps(feature(triangle, 177)))
    project = Transformer.from_crs(4326, 32644, always_xy=True).transform
    contours, spots, _, counts, _ = prepare_constraints(
        [path], box(79.25, 19.75, 79.5, 20), project
    )
    assert not contours
    assert counts["triangular_spot_heights"] == 1
    assert spots[0][1]["spot_classification"] == "triangular_symbol"
    assert spots[0][1]["symbol_type"] == "triangle"
    assert spots[0][1]["classification_method"] == "geometry_heuristic"
    assert spots[0][1]["source_geometry_sha256"] == geometry_hash(triangle)
    assert spots[0][1]["VALUE"] == 177
    assert spots[0][0].distance(Polygon(transform(project, triangle)).centroid) < 1e-8


def test_triangle_classification_precedes_clipping_in_multipart_feature(tmp_path):
    triangle = [(100, 100), (165, 100), (132.5, 165), (100, 100)]
    path = tmp_path / "source.geojsonl"
    path.write_text(
        json.dumps(feature(MultiLineString([triangle, [(0, 0), (200, 200)]]), 180))
    )
    contours, spots, _, counts, _ = prepare_constraints(
        [path], box(100, 100, 145, 160), lambda x, y: (x, y)
    )
    assert len(contours) == 1
    assert len(spots) == 1
    assert spots[0][1]["VALUE"] == 180  # Also detect heights on a contour interval.
    assert spots[0][1]["source_part_index"] == 0
    assert counts["triangular_spot_heights"] == 1


def test_reviewed_triangle_becomes_point_and_checks_source(tmp_path):
    triangle = LineString([(100, 100), (165, 100), (132.5, 165), (100, 100)])
    record = feature(triangle, 177)
    record["properties"]["id"] = "triangle"
    path = tmp_path / "triangle.geojsonl"
    path.write_text(json.dumps(record))
    review = {
        "geometry_sha256": geometry_hash(triangle),
        "source_id": "triangle",
        "value_m": 177,
        "reason": "Printed triangular height symbol verified on map",
    }
    _, unreviewed, _, _, _ = prepare_constraints(
        [path], box(0, 0, 200, 200), lambda x, y: (x, y), triangle_side_m=0
    )
    assert not unreviewed
    contours, spots, _, counts, _ = prepare_constraints(
        [path], box(0, 0, 200, 200), lambda x, y: (x, y), reviewed_spots=[review]
    )
    assert not contours
    assert len(spots) == 1
    assert spots[0][0].x == pytest.approx(132.5)
    assert spots[0][1]["VALUE"] == 177
    assert counts["reviewed_spot_heights"] == 1
    with pytest.raises(ValueError, match="no longer matches"):
        prepare_constraints(
            [path],
            box(0, 0, 200, 200),
            lambda x, y: (x, y),
            reviewed_spots=[dict(review, value_m=178)],
        )
    with pytest.raises(ValueError, match="missing from input"):
        prepare_constraints(
            [path],
            box(0, 0, 200, 200),
            lambda x, y: (x, y),
            reviewed_spots=[dict(review, geometry_sha256="not-present")],
        )


def test_conflicting_samples_rejected_without_discarding_zero(tmp_path):
    xy = np.array([[0.0, 0.0], [0.0, 0.0], [10.0, 0.0]])
    _, heights = unique_samples(xy, np.array([0.0, 0.0, 5.0]))
    assert heights.tolist() == [0.0, 5.0]
    with pytest.raises(ValueError, match="conflicting"):
        unique_samples(
            xy,
            np.array([0.0, 10.0, 5.0]),
            conflict_path=tmp_path / "conflicts.geojsonl",
            unproject=lambda x, y: (x, y),
        )
    conflict = json.loads((tmp_path / "conflicts.geojsonl").read_text())
    assert conflict["geometry"]["coordinates"] == [0, 0]
    assert conflict["properties"] == {"elevations_m": [0, 10], "range_m": 10}


def test_missing_value_is_not_sea_level(tmp_path):
    path = tmp_path / "bad.geojsonl"
    path.write_text(
        json.dumps(
            {
                "type": "Feature",
                "properties": {},
                "geometry": mapping(LineString([(0, 0), (1, 1)])),
            }
        )
    )
    with pytest.raises(ValueError, match="Invalid or missing VALUE"):
        prepare_constraints([path], box(0, 0, 1, 1), lambda x, y: (x, y))


def test_zero_length_source_line_is_counted_without_becoming_a_symbol(tmp_path):
    path = tmp_path / "source.geojsonl"
    path.write_text(json.dumps(feature(LineString([(0.5, 0.5), (0.5, 0.5)]), 120)))
    contours, spots, _, counts, _ = prepare_constraints(
        [path], box(0, 0, 1, 1), lambda x, y: (x, y)
    )
    assert contours == spots == []
    assert counts["zero_length_line_parts"] == 1


@pytest.mark.parametrize(
    "variation,expected", [("exact", 1), ("height_mismatch", 0), ("gap", 0)]
)
def test_reconstruct_split_circle_preserves_sources_and_requires_exact_join(
    tmp_path, variation, expected
):
    ring = list(Point(100, 100).buffer(9, quad_segs=16).exterior.coords)
    a, b = ring[:33], ring[32:]
    if variation == "gap":
        b[0] = (b[0][0] + 0.001, b[0][1])
    records = [
        feature(LineString(a), 529),
        feature(LineString(b), 530 if variation == "height_mismatch" else 529),
    ]
    for i, record in enumerate(records):
        record["properties"].update(id=f"arc-{i}", OBJECT_ID=f"original-{i}")
    # A different-height contour crossing the circle must stay a contour.
    records.append(feature(LineString([(0, 100), (200, 100)]), 540))
    source = tmp_path / "source.geojsonl"
    source.write_text("\n".join(json.dumps(r) for r in records))
    contours, spots, _, counts, _ = prepare_constraints(
        [source], box(0, 0, 200, 200), lambda x, y: (x, y)
    )
    assert len(spots) == expected
    if expected:
        assert len(contours) == 1
        point, properties = spots[0]
        assert point.distance(Point(100, 100)) < 1e-8
        assert (
            properties["classification_method"] == "geometry_heuristic_joined_fragments"
        )
        assert properties["elevation_m"] == 529
        assert [p["source_id"] for p in properties["source_fragments"]] == [
            "arc-0",
            "arc-1",
        ]
        assert (
            properties["source_fragments"][1]["properties"]["OBJECT_ID"] == "original-1"
        )
        assert counts["reconstructed_symbol_fragments"] == 2
    else:
        assert len(contours) >= 3


def test_reconstruct_triangle_fragments_keeps_unrelated_multipart_contour(tmp_path):
    source = tmp_path / "source.geojsonl"
    a = feature(
        MultiLineString(
            [
                [(100, 100), (165, 100)],
                [(0, 0), (200, 0)],
            ]
        ),
        177,
    )
    b = feature(LineString([(165, 100), (132.5, 165), (100, 100)]), 177)
    a["properties"]["id"], b["properties"]["id"] = "a", "b"
    source.write_text("\n".join(json.dumps(f) for f in [a, b]))
    contours, spots, _, _, _ = prepare_constraints(
        [source], box(0, 0, 200, 200), lambda x, y: (x, y)
    )
    assert len(contours) == len(spots) == 1
    assert contours[0][1]["id"] == "a"
    assert spots[0][1]["symbol_type"] == "triangle"
    assert [p["source_id"] for p in spots[0][1]["source_fragments"]] == ["a", "b"]


@pytest.mark.parametrize("form_mode", ["auto", None])
def test_small_height_fragment_is_recovered_in_final_dem_with_original_audit(
    tmp_path, form_mode
):
    unproject = Transformer.from_crs(32644, 4326, always_xy=True).transform
    sheet = transform(unproject, box(540000, 2200000, 541000, 2201000))
    index = tmp_path / "index.geojson"
    index.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"id": "test"},
                        "geometry": mapping(sheet),
                    }
                ],
            }
        )
    )
    lines = [
        (LineString([(540050, y), (540950, y)]), z)
        for y, z in [(2200050, 840), (2200300, 860), (2200700, 880), (2200950, 900)]
    ]
    lines += [
        (
            LineString(
                [
                    (540506, 2200500),
                    (540530, 2200500),
                    (540530, 2200530),
                    (540500, 2200530),
                    (540500, 2200500),
                ]
            ),
            880,
        ),
        (LineString([(540500, 2200500), (540506, 2200500)]), 1240),
        (
            LineString(Point(540550, 2200550).buffer(9, quad_segs=16).exterior.coords),
            889,
        ),
    ]
    source = tmp_path / "contours.geojsonl"
    features = [feature(transform(unproject, g), z) for g, z in lines]
    for i, record in enumerate(features):
        record["properties"]["id"] = f"source-{i}"
    source.write_text("\n".join(json.dumps(f) for f in features))
    before = source.read_bytes()
    output = tmp_path / "recovered"
    report = build_dem([source], index, "test", output, form_line_interval=form_mode)
    assert report["height_recovery"]["changed_source_parts"] == 1
    assert report["height_recovery"]["methods"] == {
        "small_fragment_between_same_height_endpoints": 1
    }
    assert report["elevation_max_m"] <= 900
    assert report["spot_holdout_contours_only"]["count"] == 1
    original = [
        json.loads(s) for s in (output / "contours.geojsonl").read_text().splitlines()
    ]
    repaired = [
        json.loads(s)
        for s in (output / "recovered_contours.geojsonl").read_text().splitlines()
    ]
    assert (
        next(
            f["properties"]["VALUE"]
            for f in original
            if f["properties"]["id"] == "source-5"
        )
        == 1240
    )
    correction = next(
        f["properties"] for f in repaired if f["properties"]["id"] == "source-5"
    )
    assert correction["VALUE"] == 880
    assert correction["height_recovery"]["original_elevation_m"] == 1240
    assert (
        json.loads((output / "spots.geojsonl").read_text())["properties"]["VALUE"]
        == 889
    )
    assert source.read_bytes() == before
    for name, options in [
        ("disabled", {"height_recovery": "none"}),
        ("too_short", {"height_repair_max_length_m": 5}),
    ]:
        with pytest.raises(ValueError, match="conflicting elevations"):
            build_dem(
                [source],
                index,
                "test",
                tmp_path / name,
                form_line_interval=form_mode,
                **options,
            )
        assert not (tmp_path / name / "dem.tif").exists()


def test_planar_dem_and_projection_end_to_end(tmp_path):
    unproject = Transformer.from_crs(32644, 4326, always_xy=True).transform
    sheet = transform(unproject, box(540000, 2200000, 541000, 2201000))
    index = tmp_path / "index.geojson"
    index.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"id": "test"},
                        "geometry": mapping(sheet),
                    }
                ],
            }
        )
    )
    path = tmp_path / "contours.geojsonl"
    records = [
        feature(transform(unproject, LineString([(x, 2200000), (x, 2201000)])), z)
        for x, z in [(540000, 0), (540500, 50), (541000, 100)]
    ]
    circle = LineString(Point(540250, 2200300).buffer(9, quad_segs=16).exterior.coords)
    triangle = LineString(
        [
            (540717.5, 2200700),
            (540782.5, 2200700),
            (540750, 2200765),
            (540717.5, 2200700),
        ]
    )
    for symbol, elevation, source_id in [(circle, 25, "c1"), (triangle, 75, "t1")]:
        record = feature(transform(unproject, symbol), str(elevation))
        record["properties"].update(id=source_id, OBJECT_ID="original", TYPE="test")
        records.append(record)
    path.write_text("\n".join(json.dumps(f) for f in records))
    report = build_dem([path], index, "test", tmp_path / "out", resolution=30)
    assert report["crs"] == "EPSG:32644"
    assert report["coverage_fraction"] > 0.95
    assert report["symbol_type_counts"] == {"circle": 1, "triangle": 1}
    points = [
        json.loads(line)
        for line in (tmp_path / "out/spots.geojsonl").read_text().splitlines()
    ]
    assert len(points) == 2
    for point, symbol_type, source_id in zip(
        points, ["circle", "triangle"], ["c1", "t1"]
    ):
        props = point["properties"]
        assert point["geometry"]["type"] == "Point"
        assert props["symbol_type"] == symbol_type
        assert props["id"] == props["source_id"] == source_id
        assert props["source_sheet"] == "test"
        assert props["OBJECT_ID"] == "original"
        assert props["TYPE"] == "test"
        assert props["elevation_m"] == props["VALUE"]
        assert len(props["source_geometry_sha256"]) == 64
    assert points[0]["properties"]["symbol_radius_m"] == pytest.approx(9)
    assert points[1]["properties"]["symbol_side_min_m"] == pytest.approx(65)
    assert (
        json.loads((tmp_path / "out/triangle_points.geojsonl").read_text()) == points[1]
    )
    with rasterio.open(tmp_path / "out/dem.tif") as raster:
        data = raster.read(1, masked=True)
        rows, cols = np.where(~data.mask)
        xs, _ = rasterio.transform.xy(raster.transform, rows, cols)
        np.testing.assert_allclose(
            data[rows, cols], (np.asarray(xs) - 540000) / 10, atol=0.002
        )
        assert raster.res == (30, 30)
        assert raster.nodata == -9999
    with pytest.raises(ValueError, match="not empty"):
        build_dem([path], index, "test", tmp_path / "out")


def test_water_flattening_leaves_incomplete_and_variable_shores():
    model = LinearNDInterpolator(
        [[0, 0], [100, 0], [0, 100], [100, 100]], [0, 10, 0, 10]
    )
    grid = from_origin(0, 100, 10, 10)
    dem = np.tile(np.arange(10) + 0.5, (10, 1))
    flattened, decisions = flatten_water(
        dem, [box(20, 20, 40, 40)], model, grid, box(0, 0, 100, 100), 10
    )
    assert decisions[0]["status"] == "flattened"
    assert np.all(flattened[6:8, 2:4] == 3)
    for polygon, status in [
        (box(-10, 10, 20, 30), "skipped_crosses_sheet"),
        (box(10, 10, 90, 90), "skipped_variable_shoreline"),
    ]:
        flattened, decisions = flatten_water(
            dem, [polygon], model, grid, box(0, 0, 100, 100), 10
        )
        assert decisions[0]["status"] == status
        np.testing.assert_array_equal(flattened, dem)


def test_final_dem_includes_both_symbols_after_contour_only_holdout(tmp_path):
    unproject = Transformer.from_crs(32644, 4326, always_xy=True).transform
    sheet = transform(unproject, box(540000, 2200000, 541020, 2201020))
    index = tmp_path / "index.geojson"
    index.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"id": "test"},
                        "geometry": mapping(sheet),
                    }
                ],
            }
        )
    )
    records = [
        feature(transform(unproject, LineString([(x, 2200000), (x, 2201020)])), 0)
        for x in [540000, 541020]
    ]
    # Place the two symbol centres at exact output pixel centres. Both differ
    # from the flat contour-only surface, so a final model omitting them fails.
    circle_centre = Point(540315, 2200395)
    triangle_centre = Point(540795, 2200785)
    circle = LineString(circle_centre.buffer(9, quad_segs=16).exterior.coords)
    x, y = triangle_centre.x, triangle_centre.y
    triangle = LineString(
        [
            (x - 32.5, y - 65 / 3),
            (x + 32.5, y - 65 / 3),
            (x, y + 130 / 3),
            (x - 32.5, y - 65 / 3),
        ]
    )
    records.extend(
        [
            feature(transform(unproject, circle), 100),
            feature(transform(unproject, triangle), 60),
        ]
    )
    source = tmp_path / "contours.geojsonl"
    source.write_text("\n".join(json.dumps(f) for f in records))
    output = tmp_path / "dem"
    report = build_dem([source], index, "test", output)
    assert report["spot_holdout_contours_only"]["count"] == 2
    assert report["spot_holdout_contours_only"]["bias_m"] == pytest.approx(-80)
    with rasterio.open(output / "dem.tif") as raster:
        values = list(
            raster.sample([(p.x, p.y) for p in [circle_centre, triangle_centre]])
        )
        assert [v[0] for v in values] == pytest.approx([100, 60], abs=0.001)


def test_conflicting_contours_preserve_classified_points_without_a_dem(tmp_path):
    unproject = Transformer.from_crs(32644, 4326, always_xy=True).transform
    sheet = feature(transform(unproject, box(540000, 2200000, 541000, 2201000)), 0)
    sheet["properties"] = {"id": "test"}
    index = tmp_path / "index.geojson"
    index.write_text(json.dumps({"type": "FeatureCollection", "features": [sheet]}))
    shared = transform(unproject, LineString([(540100, 2200100), (540900, 2200900)]))
    circle = transform(
        unproject, LineString(Point(540500, 2200300).buffer(9).exterior.coords)
    )
    source = tmp_path / "source.geojsonl"
    source.write_text(
        "\n".join(
            json.dumps(f)
            for f in [feature(shared, 0), feature(shared, 20), feature(circle, 7)]
        )
    )
    output = tmp_path / "out"
    with pytest.raises(ValueError, match="conflicting"):
        build_dem([source], index, "test", output)
    point = json.loads((output / "spots.geojsonl").read_text())
    assert point["properties"]["symbol_type"] == "circle"
    assert point["properties"]["elevation_m"] == 7
    assert (
        json.loads((output / "classification.json").read_text())["counts"][
            "spot_heights"
        ]
        == 1
    )
    assert (output / "conflicting_constraints.geojsonl").stat().st_size > 0
    assert not (output / "dem.tif").exists()
    assert not (output / "report.json").exists()


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
