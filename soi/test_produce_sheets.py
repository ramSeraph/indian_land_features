# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "rasterio>=1.4", "scipy>=1.14", "shapely>=2.1", "gh_release_tools==0.1.14", "pytest>=8"]
# ///
import json

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from contours_to_dem import write_raster
from produce_sheets import package_sheet, retain_height_points


def test_delivery_contains_only_pixel_identical_cog(tmp_path):
    job = tmp_path / "job"
    (job / "dem").mkdir(parents=True)
    data = np.array([[100, 120], [140, np.nan]], dtype="float32")
    write_raster(
        job / "dem/dem.tif", data, from_origin(300000, 2200000, 30, 30), 32644, "test"
    )
    (job / "dem/diagnostic.json").write_text("{}")
    destination = tmp_path / "delivery"
    assets = package_sheet("56M_9", job, destination, {"status": "built"})
    assert [a["name"] for a in assets] == ["SOI_DEM_30m_56M_9.tif"]
    with rasterio.open(destination / assets[0]["name"]) as src:
        np.testing.assert_array_equal(src.read(1, masked=True).filled(np.nan), data)
        assert src.tags(ns="IMAGE_STRUCTURE")["LAYOUT"] == "COG"


def test_failed_dem_has_no_release_asset_but_keeps_height_points(tmp_path):
    job = tmp_path / "job"
    dem = job / "dem"
    dem.mkdir(parents=True)
    points = "".join(
        json.dumps(
            {
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [79.5, 19.5]},
                "properties": {"symbol_type": symbol, "elevation_m": height},
            }
        )
        + "\n"
        for symbol, height in [("circle", 177), ("triangle", 200)]
    )
    (dem / "spots.geojsonl").write_text(points)
    assert (
        package_sheet("56M_9", job, tmp_path / "delivery", {"status": "failed"}) == []
    )
    retain_height_points("56M_9", dem, tmp_path / "points")
    assert (tmp_path / "points/56M_9.geojsonl").read_text() == points
    assert not list(tmp_path.rglob("*.zip"))


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
