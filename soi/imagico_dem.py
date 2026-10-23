# SPDX-License-Identifier: GPL-3.0-or-later
# SOI raster adapter for the Python adaptation of Christoph Hormann's method.
# See LICENSE.imagico.txt and imagico_provenance.json.
"""Prepare vector constraints for Imagico interpolation and sample its grid."""

import hashlib
import math
from pathlib import Path

import numpy as np
from scipy import ndimage
from shapely import segmentize

from imagico_interpolate import interpolate


def _pixels(coordinates, shape, grid_transform):
    x, y = np.asarray(coordinates, dtype=float).T
    inverse = ~grid_transform
    cols = inverse.a * x + inverse.b * y + inverse.c
    rows = inverse.d * x + inverse.e * y + inverse.f
    height, width = shape
    inside = (
        (cols >= -1e-6)
        & (rows >= -1e-6)
        & (cols <= width + 1e-6)
        & (rows <= height + 1e-6)
    )
    col = np.floor(np.clip(cols[inside], 0, width - 1e-7)).astype(int)
    row = np.floor(np.clip(rows[inside], 0, height - 1e-7)).astype(int)
    distance2 = (cols[inside] - col - 0.5) ** 2 + (rows[inside] - row - 0.5) ** 2
    return row * width + col, distance2


def rasterize_constraints(contours, spots, unknown_contours, shape, grid_transform):
    """Use the closest sampled line per pixel; exact point heights take priority.

    Different known values can occupy one pixel at finite resolution. This is
    reported separately from contradictory vector junctions, which the caller must
    check first. Different point heights in one pixel are refused, never averaged.
    """
    if (
        grid_transform.a <= 0
        or grid_transform.b
        or grid_transform.d
        or not math.isclose(grid_transform.a, -grid_transform.e)
    ):
        raise ValueError("Imagico interpolation requires a north-up square grid")
    resolution = grid_transform.a
    indices, distances, elevations = [], [], []
    for line, props in contours:
        coords = np.asarray(segmentize(line, resolution / 4).coords)[:, :2]
        cells, delta = _pixels(coords, shape, grid_transform)
        indices.extend(cells)
        distances.extend(delta)
        elevations.extend([float(props["VALUE"])] * len(cells))
    cells = np.asarray(indices, dtype=int)
    z = np.asarray(elevations)
    known = np.full(shape, np.nan)
    low = np.full(known.size, np.inf)
    high = np.full(known.size, -np.inf)
    if len(cells):
        np.minimum.at(low, cells, z)
        np.maximum.at(high, cells, z)
        # Elevation breaks exact geometric ties deterministically, independently
        # of the order in which release features were read.
        order = np.lexsort((z, np.asarray(distances), cells))
        _, first = np.unique(cells[order], return_index=True)
        selected = order[first]
        known.ravel()[cells[selected]] = z[selected]
    collisions = np.isfinite(low) & (high - low > 1e-7)
    point_values = {}
    point_contour_collisions = 0
    for point, props in spots:
        point_index, _ = _pixels([[point.x, point.y]], shape, grid_transform)
        if not len(point_index):
            raise ValueError("A height point lies outside the interpolation grid")
        cell, height = int(point_index[0]), float(props["VALUE"])
        if cell in point_values and not math.isclose(
            point_values[cell], height, abs_tol=1e-7, rel_tol=0
        ):
            raise ValueError(
                "Different point heights occupy one grid cell; use a finer resolution"
            )
        if (
            np.isfinite(known.ravel()[cell])
            and abs(known.ravel()[cell] - height) > 1e-7
        ):
            point_contour_collisions += 1
        point_values[cell] = height
    for cell, height in point_values.items():
        known.ravel()[cell] = height
    unknown = np.zeros(shape, dtype=np.int32)
    group_pixels = {}
    for line, props in unknown_contours:
        group = int(props["unknown_height_id"])
        if group <= 0:
            raise ValueError("Unknown contour IDs must be positive")
        coords = np.asarray(segmentize(line, resolution / 4).coords)[:, :2]
        selected, _ = _pixels(coords, shape, grid_transform)
        group_pixels.setdefault(group, set()).update(int(i) for i in selected)
    for group, selected in sorted(group_pixels.items()):
        if not selected:
            raise ValueError(f"Unknown contour {group} has no grid pixels")
        indices = np.fromiter(selected, dtype=int)
        if (unknown.ravel()[indices] != 0).any():
            raise ValueError(
                "Different unknown contours occupy one grid cell; use a finer resolution"
            )
        unknown.ravel()[indices] = group
    report = {
        "line_sampling_spacing_m": resolution / 4,
        "contour_pixel_selection": "nearest_sample_to_pixel_centre",
        "contour_collision_pixels": int(collisions.sum()),
        "maximum_contour_pixel_height_range_m": float(np.max((high - low)[collisions]))
        if collisions.any()
        else 0,
        "height_points_used": len(spots),
        "distinct_point_pixels": len(point_values),
        "point_contour_collision_count": point_contour_collisions,
        "point_pixel_priority": True,
        "unknown_contour_groups": len(group_pixels),
    }
    return known, unknown, report


def sample_grid(grid, grid_transform, coordinates):
    coordinates = np.asarray(coordinates)
    inverse = ~grid_transform
    x, y = coordinates[:, 0], coordinates[:, 1]
    cols = inverse.a * x + inverse.b * y + inverse.c
    rows = inverse.d * x + inverse.e * y + inverse.f
    return ndimage.map_coordinates(
        grid, [rows - 0.5, cols - 0.5], order=1, mode="constant", cval=np.nan
    )


def build_grid(
    contours,
    spots,
    unknown_contours,
    shape,
    grid_transform,
    *,
    progress=None,
    **options,
):
    known, unknown, raster_report = rasterize_constraints(
        contours, spots, unknown_contours, shape, grid_transform
    )
    result = interpolate(known, unknown, progress=progress, **options)
    result.report["rasterization"] = raster_report
    result.report["implementation"] = {
        "license": "GPL-3.0-or-later",
        "provenance": "imagico_provenance.json",
        "files": [
            {
                "path": name,
                "sha256": hashlib.sha256(
                    Path(__file__).with_name(name).read_bytes()
                ).hexdigest(),
            }
            for name in (
                "imagico_interpolate.py",
                "imagico_dem.py",
                "imagico_provenance.json",
                "LICENSE.imagico.txt",
            )
        ],
    }
    return result, known, unknown
