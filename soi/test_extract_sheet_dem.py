# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "rasterio>=1.4", "scipy>=1.14", "shapely>=2.1", "pytest>=8"]
# ///
import json
import subprocess

import pytest
from pyproj import Transformer
from shapely.geometry import LineString, box, mapping

import extract_sheet_dem as pipeline
from contours_to_dem import transform


def test_sheet_pipeline_and_cached_provenance_without_toposheet(tmp_path, monkeypatch):
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
    contours = [
        {
            "type": "Feature",
            "properties": {"VALUE": z},
            "geometry": mapping(
                transform(unproject, LineString([(x, 2200000), (x, 2201000)]))
            ),
        }
        for x, z in [(540000, 0), (540500, 50), (541000, 100)]
    ]
    calls = []

    def iomaps(command, check):
        assert check
        assert command[command.index("--http-client") + 1] == "curl"
        assert command[command.index("-f") + 1] == str(index)
        assert command[command.index("--pick-filter-feature-kv") + 1] == "id=test"
        source = command[command.index("-s") + 1]
        calls.append(source)
        target = pipeline.Path(command[command.index("-o") + 1])
        records = contours if source == pipeline.CONTOURS else []
        target.write_text("\n".join(json.dumps(f) for f in records))

    monkeypatch.setattr(pipeline.shutil, "which", lambda executable: "/fake/iomaps")
    monkeypatch.setattr(pipeline.subprocess, "run", iomaps)
    cache = tmp_path / "cache"
    output = tmp_path / "dem"
    report = pipeline.run_sheet(
        "test", index, output, cache, http_client="curl",
        height_repair_max_length_m=7.0, height_repair_max_fraction=0.05,
    )
    assert calls == [pipeline.CONTOURS, pipeline.WATER]
    assert report["coverage_fraction"] > 0.95
    assert report["height_recovery"]["fragment_max_length_m"] == 7.0
    assert report["height_recovery"]["fragment_max_fraction"] == 0.05
    assert (output / "dem_water_flattened.tif").exists()
    assert (output / "triangle_points.geojsonl").read_text() == ""
    extraction = json.loads((output / "extraction.json").read_text())
    assert [r["features"] for r in extraction["layers"]] == [3, 0]
    assert all(r["origin"] == "iomaps" for r in extraction["layers"])
    # Cache reuse needs neither the iomaps binary nor any toposheet.
    monkeypatch.setattr(pipeline.shutil, "which", lambda executable: None)
    pipeline.run_sheet("test", index, tmp_path / "reused", cache, reuse=True)
    assert len(calls) == 2
    cached = cache / "test/contours.geojsonl"
    cached.write_text(cached.read_text() + "\n")
    with pytest.raises(ValueError, match="does not match its metadata"):
        pipeline.run_sheet("test", index, tmp_path / "changed", cache, reuse=True)


@pytest.mark.parametrize("failure", ["interrupted", "empty"])
def test_failed_extraction_does_not_publish_reusable_file(
    tmp_path, monkeypatch, failure
):
    def iomaps(command, check):
        target = pipeline.Path(command[command.index("-o") + 1])
        target.write_text("incomplete" if failure == "interrupted" else "")
        if failure == "interrupted":
            raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr(pipeline.shutil, "which", lambda executable: "/fake/iomaps")
    monkeypatch.setattr(pipeline.subprocess, "run", iomaps)
    target = tmp_path / "contours.geojsonl"
    with pytest.raises((subprocess.CalledProcessError, ValueError)):
        pipeline.extract_layer(
            "iomaps", pipeline.CONTOURS, "test", "hash", tmp_path / "index", target
        )
    assert not target.exists()
    assert not target.with_suffix(".metadata.json").exists()


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
