# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "rasterio>=1.4", "scipy>=1.14", "shapely>=2.1", "pytest>=8"]
# ///
# SPDX-License-Identifier: GPL-3.0-or-later
import json

import numpy as np
import pytest
import rasterio
from pyproj import Transformer
from rasterio.transform import from_origin
from shapely.geometry import LineString, Point, box, mapping

from contours_to_dem import build_dem, transform
from imagico_dem import rasterize_constraints


def inputs(tmp_path, lines):
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
    source = tmp_path / "source.geojsonl"
    source.write_text(
        "\n".join(
            json.dumps(
                {
                    "type": "Feature",
                    "properties": {"id": f"source-{i}", "VALUE": z},
                    "geometry": mapping(transform(unproject, g)),
                }
            )
            for i, (g, z) in enumerate(lines)
        )
    )
    return source, index


def test_default_engine_holds_out_then_restores_circle_and_triangle_points(tmp_path):
    lines = [
        (LineString([(x, 2200000), (x, 2201000)]), 0) for x in [540000, 540500, 541000]
    ]
    circle = LineString(Point(540315, 2200395).buffer(9, quad_segs=16).exterior.coords)
    x, y = 540795, 2200785
    triangle = LineString(
        [
            (x - 32.5, y - 65 / 3),
            (x + 32.5, y - 65 / 3),
            (x, y + 130 / 3),
            (x - 32.5, y - 65 / 3),
        ]
    )
    source, index = inputs(tmp_path, lines + [(circle, 100), (triangle, 60)])
    output = tmp_path / "dem"
    report = build_dem([source], index, "test", output)
    assert report["method"] == "imagico_python"
    assert report["symbol_type_counts"] == {"circle": 1, "triangle": 1}
    assert report["spot_holdout_contours_only"]["count"] == 2
    assert report["spot_holdout_contours_only"]["bias_m"] == pytest.approx(-80)
    assert report["interpolation_engine"]["rasterization"]["height_points_used"] == 2
    validation = json.loads((output / "imagico_validation_report.json").read_text())
    assert validation["rasterization"]["height_points_used"] == 0
    assert report["interpolation_engine"]["known_constraints_preserved"]
    with rasterio.open(output / "dem.tif") as raster:
        assert raster.tags()["method"] == "imagico_python"
        heights = [v[0] for v in raster.sample([(540315, 2200395), (540795, 2200785)])]
        assert heights == pytest.approx([100, 60], abs=1e-3)
        validated_final = raster.read(1)
    production = tmp_path / "production"
    final_report = build_dem([source], index, "test", production, validate_points=False)
    assert not (production / "imagico_validation_report.json").exists()
    assert not (production / "dem_contours_only.tif").exists()
    assert final_report["interpolation_engine"]["rasterization"]["height_points_used"] == 2
    with rasterio.open(production / "dem.tif") as raster:
        np.testing.assert_array_equal(raster.read(1), validated_final)


def test_suspect_height_becomes_unknown_instead_of_a_fixed_replacement(tmp_path):
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
        # This point occupies the unknown fragment's pixel. Its height must
        # affect the final estimate even though it was absent in validation.
        (Point(540505, 2200499), 889),
    ]
    source, index = inputs(tmp_path, lines)
    original = source.read_bytes()
    output = tmp_path / "dem"
    report = build_dem([source], index, "test", output)
    unknown = [
        json.loads(s)
        for s in (output / "unknown_height_contours.geojsonl").read_text().splitlines()
    ]
    assert len(unknown) == 1
    props = unknown[0]["properties"]
    assert props["VALUE"] is None
    assert props["original_elevation_m"] == 1240
    assert report["height_recovery"]["application"] == "unassign_then_interpolate"
    estimates = report["interpolation_engine"]["unknown_height_estimates"]
    assert len(estimates) == 1
    assert 840 <= estimates[0]["estimated_elevation_m"] <= 900
    validation = json.loads((output / "imagico_validation_report.json").read_text())
    assert validation["rasterization"]["height_points_used"] == 0
    assert report["interpolation_engine"]["rasterization"]["height_points_used"] == 1
    assert estimates[0]["estimated_elevation_m"] != pytest.approx(
        validation["unknown_height_estimates"][0]["estimated_elevation_m"], abs=1e-5
    )
    grids = np.load(output / "interpolation_inputs.npz")
    assert np.nanmax(grids["known_heights"]) <= 900
    assert source.read_bytes() == original
    recovered = [
        json.loads(s)
        for s in (output / "recovered_contours.geojsonl").read_text().splitlines()
    ]
    assigned = [f for f in recovered if f["properties"].get("unknown_height_id")]
    assert assigned[0]["properties"]["VALUE"] == estimates[0]["estimated_elevation_m"]


def test_raster_collisions_are_deterministic_and_keep_point_heights():
    lines = [
        (LineString([(1, 1), (1, 19)]), {"VALUE": 10}),
        (LineString([(2, 1), (2, 19)]), {"VALUE": 20}),
    ]
    points = [(Point(2, 2), {"VALUE": 15})]
    a, _, report = rasterize_constraints(
        lines, points, [], (5, 5), from_origin(0, 50, 10, 10)
    )
    b, _, _ = rasterize_constraints(
        lines[::-1], points, [], (5, 5), from_origin(0, 50, 10, 10)
    )
    assert np.array_equal(a, b, equal_nan=True)
    assert a[4, 0] == 15
    assert report["contour_collision_pixels"] > 0
    with pytest.raises(ValueError, match="Different point heights"):
        rasterize_constraints(
            lines,
            points + [(Point(3, 3), {"VALUE": 16})],
            [],
            (5, 5),
            from_origin(0, 50, 10, 10),
        )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
