# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "rasterio>=1.4", "scipy>=1.14", "pyproj>=3.6", "shapely>=2.1", "pytest>=8"]
# ///
import json
import math
import shutil

import numpy as np
import pytest
import rasterio
from pyproj import Transformer
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from shapely import segmentize
from shapely.geometry import box, mapping
from shapely.ops import unary_union

from contours_to_dem import transform, write_raster
from overlap_cache import OverlapProduction
from sheet_overlaps import SheetIndex, correct_sheet, create_patch


@pytest.fixture
def four_sheets(tmp_path):
    polygons = {
        "nw": box(79.495, 19.75, 79.5, 19.755),
        "ne": box(79.5, 19.75, 79.505, 19.755),
        "sw": box(79.495, 19.745, 79.5, 19.75),
        "se": box(79.5, 19.745, 79.505, 19.75),
    }
    index_path = tmp_path / "index.json"
    index_path.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {"type": "Feature", "properties": {"id": s}, "geometry": mapping(g)}
                    for s, g in polygons.items()
                ],
            }
        )
    )
    index = SheetIndex(index_path, polygons)
    project = Transformer.from_crs(4326, 32644, always_xy=True).transform
    region = transform(project, segmentize(unary_union(list(polygons.values())), 0.002))
    west, south, east, north = region.bounds
    west, south = math.floor(west / 30) * 30, math.floor(south / 30) * 30
    east, north = math.ceil(east / 30) * 30, math.ceil(north / 30) * 30
    affine = from_origin(west, north, 30, 30)
    grid_shape = (round((north - south) / 30), round((east - west) / 30))
    folders = {}
    for s, offset in zip(polygons, [-4, 2, -2, 4], strict=True):
        folder = tmp_path / s
        folder.mkdir()
        inside = geometry_mask(
            [mapping(transform(project, segmentize(polygons[s], 0.002)))],
            grid_shape,
            affine,
            invert=True,
        )
        data = np.where(inside, 100.0 + offset, np.nan)
        known = np.full(grid_shape, np.nan)
        rr, cc = np.nonzero(inside)
        row, col = round(rr.mean()), round(cc.mean())
        known[row, col] = data[row, col] = 100
        write_raster(folder / "dem.tif", data, affine, 32644, "test")
        np.savez_compressed(
            folder / "interpolation_inputs.npz",
            known_heights=known,
            unknown_ids=np.zeros(grid_shape, dtype=np.int32),
        )
        (folder / "report.json").write_text(
            json.dumps({"sheet": s, "source_files": []})
        )
        folders[s] = folder
    return index_path, index, folders


def test_canonical_jobs_and_diagonal_only_neighbours(four_sheets):
    path, index, folders = four_sheets
    assert sum(j["kind"] == "edge" for j in index.jobs.values()) == 4
    assert sum(j["kind"] == "corner" for j in index.jobs.values()) == 1
    assert index.neighbourhood("nw") == sorted(folders)
    assert not SheetIndex(path, ["nw", "se"]).jobs
    corner = next(j for j in index.jobs.values() if j["kind"] == "corner")
    assert len(index.incident_edges(corner)) == 4


def test_corner_uses_shared_edge_patches_and_preserves_fixed_cells(
    four_sheets, tmp_path
):
    _, index, folders = four_sheets
    patches = {}
    for key, job in index.jobs.items():
        if job["kind"] == "edge":
            path = tmp_path / key
            create_patch(
                job, {s: folders[s] for s in job["members"]}, index, path, width_m=120
            )
            patches[key] = path
    key, corner = next((k, j) for k, j in index.jobs.items() if j["kind"] == "corner")
    path = tmp_path / key
    report = create_patch(
        corner,
        folders,
        index,
        path,
        prior_patches=[patches[k] for k in index.incident_edges(corner)],
        width_m=120,
    )
    patches[key] = path
    assert len(report["prior_edge_patches"]) == 4
    for s in folders:
        own = {k: p for k, p in patches.items() if s in index.jobs[k]["members"]}
        statuses = {k: {"status": "applied"} for k in own}
        a, b = tmp_path / f"{s}-a", tmp_path / f"{s}-b"
        shutil.copytree(folders[s], a)
        shutil.copytree(folders[s], b)
        result = correct_sheet(s, a, index, own, statuses, width_m=120)
        correct_sheet(
            s, b, index, dict(reversed(list(own.items()))), statuses, width_m=120
        )
        with rasterio.open(a / "dem.tif") as x, rasterio.open(b / "dem.tif") as y:
            np.testing.assert_array_equal(x.read(1), y.read(1))
        assert result["seam_reconciled"]
        assert result["fixed_constraint_pixels_checked"] == 1
        assert result["outside_patches_unchanged"]


def test_cross_utm_job_is_rejected_without_creating_output(four_sheets, tmp_path):
    _, index, folders = four_sheets
    job = next(j for j in index.jobs.values() if j["kind"] == "edge")
    path = folders[job["members"][1]] / "dem.tif"
    with rasterio.open(path, "r+") as src:
        src.crs = "EPSG:32643"
    target = tmp_path / "rejected"
    with pytest.raises(ValueError, match="Cross-UTM"):
        create_patch(
            job, {s: folders[s] for s in job["members"]}, index, target, width_m=120
        )
    assert not target.exists()


def test_production_cache_handles_missing_neighbour_and_resume(four_sheets, tmp_path):
    index_path, _, folders = four_sheets
    calls = []

    def build(sheet, job):
        calls.append(sheet)
        (job / "build.log").write_text("test build\n")
        if sheet == "se":
            return {
                "status": "failed",
                "sheet": sheet,
                "error": "conflicting source heights",
            }
        shutil.copytree(folders[sheet], job / "dem")
        return {"status": "built", "sheet": sheet}

    root = tmp_path / "production"
    root.mkdir()
    runner = OverlapProduction(
        root, index_path, folders, "test-run", build, width_m=120
    )
    result = runner.prepare("nw", root / "job-a")
    assert result["overlap"]["status"] == "partial"
    assert not result["overlap"]["seam_reconciled"]
    assert result["overlap"]["neighbour_failures"]["se"]["status"] == "failed"
    # Repeated ownership of the same corner uses the cached canonical solution.
    repeat = runner.prepare("ne", root / "job-b")
    assert len(calls) == 4
    assert repeat["overlap"]["status"] == "partial"
    assert not (root / "job-a/overlaps").exists()
    assert list((root / "patch-cache").glob("corner-*/inputs.npz"))


def test_low_disk_space_discards_only_unpinned_reproducible_cache(
    four_sheets, tmp_path, monkeypatch
):
    from types import SimpleNamespace

    index_path, _, folders = four_sheets
    root = tmp_path / "cache-pressure"
    root.mkdir()
    runner = OverlapProduction(root, index_path, folders, "cache-test", lambda *_: None)
    for name in ("nw", "ne"):
        directory = runner.sheets / name
        directory.mkdir()
        (directory / "cache_record.json").write_text("{}")
    runner.active_sheets = {"nw"}
    monkeypatch.setattr(
        "overlap_cache.shutil.disk_usage",
        lambda _: SimpleNamespace(
            free=(1 if (runner.sheets / "ne").exists() else 4) * 1024**3
        ),
    )
    runner.space_check()
    assert (runner.sheets / "nw").exists()
    assert not (runner.sheets / "ne").exists()
    runner.trim()
    assert not runner.active_sheets


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
