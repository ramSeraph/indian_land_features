"""Sanity checks on exported SOI DEMs, without a toposheet raster."""

import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.features import geometry_mask
from shapely import segmentize
from shapely.geometry import mapping, shape

from contours_to_dem import geometry_hash, get_sheet, read_features, transform


def file_hash(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def source_signature(path):
    result = Counter()
    for feature in read_features(path):
        props = dict(feature["properties"])
        if isinstance(props.get("bbox"), str):
            props["bbox"] = json.loads(props["bbox"])
        result[
            (
                geometry_hash(shape(feature["geometry"])),
                json.dumps(props, sort_keys=True),
            )
        ] += 1
    return result


def validate_sheet(directory, index, cached_extract=None):
    if not __debug__:
        raise RuntimeError("Run validation without Python optimization flags")
    directory = Path(directory)
    report = json.loads((directory / "report.json").read_text())
    source = Path(report["source_files"][0]["path"])
    assert file_hash(source) == report["source_files"][0]["sha256"], (
        "Source changed after interpolation"
    )
    result = {"sheet": report["sheet"], "dem_sha256": file_hash(directory / "dem.tif")}
    if cached_extract is not None and Path(cached_extract).exists():
        current, previous = source_signature(source), source_signature(cached_extract)
        assert current == previous, (
            "Local extraction differs from the previous remote extract"
        )
        result["extraction_parity"] = {
            "status": "identical",
            "rows": current.total(),
            "previous": str(cached_extract),
        }
    else:
        result["extraction_parity"] = {"status": "no_previous_extract"}
    with rasterio.open(directory / "dem.tif") as src:
        assert src.count == 1 and src.dtypes == ("float32",)
        assert src.crs.to_string() == report["crs"]
        assert src.res == (report["resolution_m"], report["resolution_m"])
        data = src.read(1, masked=True).filled(np.nan)
        assert not np.isinf(data).any()
        project = Transformer.from_crs(4326, src.crs, always_xy=True).transform
        sheet = transform(project, segmentize(get_sheet(index, report["sheet"]), 0.002))
        inside = geometry_mask([mapping(sheet)], src.shape, src.transform, invert=True)
        finite = np.isfinite(data)
        assert not finite[~inside].any(), "DEM leaks outside its sheet"
        assert int(finite.sum()) == report["valid_cells"] > 0
        assert int(inside.sum()) == report["sheet_cells"]
        assert np.isclose(finite.sum() / inside.sum(), report["coverage_fraction"])
        with np.load(directory / "interpolation_inputs.npz") as inputs:
            known = inputs["known_heights"]
            assert known.shape == src.shape
            fixed = np.isfinite(known) & inside
            assert np.array_equal(data[fixed], known[fixed].astype("float32")), (
                "Known constraint pixels were altered"
            )
            result["fixed_constraint_pixels_checked"] = int(fixed.sum())
        points = list(read_features(directory / "spots.geojsonl"))
        errors = []
        for feature in points:
            x, y = project(*feature["geometry"]["coordinates"][:2])
            col, row = ~src.transform * (x, y)
            assert (
                -1e-6 <= col <= src.width + 1e-6 and -1e-6 <= row <= src.height + 1e-6
            )
            col, row = (
                int(np.clip(col, 0, src.width - 1)),
                int(np.clip(row, 0, src.height - 1)),
            )
            expected = np.float32(feature["properties"]["elevation_m"])
            assert np.isfinite(data[row, col]), (
                "A height point is masked in the final DEM"
            )
            assert data[row, col] == expected, (
                "A height point was not restored in the final DEM"
            )
            errors.append(
                float(data[row, col]) - float(feature["properties"]["elevation_m"])
            )
        rasterization = report["interpolation_engine"]["rasterization"]
        assert rasterization["height_points_used"] == len(points)
        symbols = Counter(f["properties"]["symbol_type"] for f in points)
        assert dict(symbols) == report["symbol_type_counts"]
        result.update(
            status="passed",
            height_points_checked=len(points),
            symbol_counts=dict(symbols),
            max_point_float32_error_m=max(map(abs, errors), default=0),
            coverage_fraction=report["coverage_fraction"],
            contour_interval=report["form_line_detection"]["sheet_interval"],
            contour_pixel_collisions=rasterization["contour_collision_pixels"],
            elevation_min_m=float(np.nanmin(data)),
            elevation_max_m=float(np.nanmax(data)),
            seam_reconciled=report["seam_reconciled"],
        )
    (directory / "sanity_checks.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    return result
