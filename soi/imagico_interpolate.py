# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2006-2007 Christoph Hormann
# Python adaptation, 2026, for the indian_land_features project.
# Derived from contour_interpolate 0.3, DemProc.cpp. See imagico_provenance.json
# and LICENSE.imagico.txt for provenance, modifications and license terms.
"""Contour interpolation following Imagico's stages, using NumPy and SciPy.

Coordinates are square grid pixels. NaN means no constraint; zero and negative
heights are valid. This is an adaptation, not a bit-for-bit CImg reproduction.
In particular, explicit masks and exact restoration protect measured heights.
"""

from dataclasses import dataclass
import math
import time

import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve


@dataclass
class InterpolationResult:
    elevation: np.ndarray
    constraints: np.ndarray
    unknown_heights: list
    report: dict


def _slices(length, offset):
    return (
        (slice(0, length - offset), slice(offset, length))
        if offset >= 0
        else (slice(-offset, length), slice(0, length + offset))
    )


def directional_fill(fixed, radius=120):
    """First known pixel in each of 32 directions, inverse-distance weighted.

    Unreached cells stay NaN; they are never confused with a valid sea-level value.
    Integer ray offsets follow DirSearchFill in the original implementation.
    """
    fixed = np.asarray(fixed, dtype=float)
    known = np.isfinite(fixed)
    height, width = fixed.shape
    weighted = np.zeros_like(fixed)
    weights = np.zeros_like(fixed)
    radius = min(radius, math.ceil(math.hypot(width, height)) + 1)
    for direction in range(32):
        remaining = ~known
        angle = 2 * math.pi * direction / 32
        previous = (0, 0)
        for step in range(1, radius + 1):
            dx, dy = int(math.sin(angle) * step), int(math.cos(angle) * step)
            if (dx, dy) == previous:
                continue
            previous = dx, dy
            if abs(dx) >= width or abs(dy) >= height:
                break
            dst_y, src_y = _slices(height, dy)
            dst_x, src_x = _slices(width, dx)
            active = remaining[dst_y, dst_x]
            hits = active & known[src_y, src_x]
            weight = 1 / math.hypot(dx, dy)
            weighted[dst_y, dst_x][hits] += fixed[src_y, src_x][hits] * weight
            weights[dst_y, dst_x][hits] += weight
            active[hits] = False
            # Pixels outside this rectangle have left the image along this ray.
            if dy > 0:
                remaining[height - dy :] = False
            elif dy < 0:
                remaining[:-dy] = False
            if dx > 0:
                remaining[:, width - dx :] = False
            elif dx < 0:
                remaining[:, :-dx] = False
            if not remaining.any():
                break
    result = np.full_like(fixed, np.nan)
    np.divide(weighted, weights, out=result, where=weights > 0)
    result[known] = fixed[known]
    return result


def _padding_surface(surface):
    """Nearest values are numerical padding only, not newly supported output."""
    valid = np.isfinite(surface)
    if not valid.any():
        raise ValueError("No supported interpolation cells")
    if valid.all():
        return surface.copy()
    nearest = ndimage.distance_transform_edt(
        ~valid, return_distances=False, return_indices=True
    )
    return surface[tuple(nearest)]


def _laplace_correction(fixed, sigma, iterations, radius, max_error):
    """Weighted spreading followed by repeated blur/constraint restoration."""
    valid = np.isfinite(fixed)
    result = fixed.copy()
    yy, xx = np.mgrid[-radius : radius + 1, -radius : radius + 1]
    distance = np.hypot(xx, yy)
    kernel = 1 / np.maximum(distance, 1)
    kernel[radius, radius] = 5
    for _ in range(5):
        available = np.isfinite(result)
        sums = fftconvolve(np.where(available, result, 0), kernel, mode="same")
        weights = fftconvolve(available.astype(float), kernel, mode="same")
        can_fill = ~valid & (weights > 1e-8)
        result[can_fill] = sums[can_fill] / weights[can_fill]
        result[valid] = fixed[valid]
    # An unpropagated correction means no correction, not sea-level terrain.
    result = np.nan_to_num(result, nan=0)
    for _ in range(iterations):
        result[valid] = fixed[valid]
        result = ndimage.gaussian_filter(result, sigma, mode="nearest")
    error = float(np.max(np.abs(result[valid] - fixed[valid]))) if valid.any() else 0
    result[valid] = fixed[valid]
    return result, max(1.0, 0.75 * error / max_error)


def _refine(initial, fixed, iterations, max_error, final=False):
    valid = np.isfinite(fixed)
    surface = ndimage.gaussian_filter(initial, 4.8, mode="nearest")
    sigma = 4.2 if final else 2.2
    residual = np.where(valid, fixed - surface, np.nan)
    correction, scale = _laplace_correction(
        residual, sigma, iterations, 4 if final else 5, max_error
    )
    surface += ndimage.gaussian_filter(correction, sigma, mode="nearest") * (
        0.95 if final else 0.85
    )
    snapshot = surface.copy()
    scale = 1.0
    for step in range(10 if final else 2):
        sigma = (3.8 - 3.2 * step / 10) / 2 / scale if final else 1.5 - 0.8 * step
        residual = np.where(valid, fixed - surface, np.nan)
        correction, scale = _laplace_correction(
            residual, sigma, iterations, 4 if final else 15, max_error
        )
        surface += (
            ndimage.gaussian_filter(correction, sigma / scale, mode="nearest") * 0.85
        )
    return surface, snapshot


_NEIGHBOURS = ((0, 1), (0, -1), (1, 0), (-1, 0), (-1, -1), (1, 1), (-1, 1), (1, -1))


def _extrema_constraints(surface, fixed, snapshot, support, step=40):
    """Trace local extrema and estimate summit/depression constraints.

    Known heights at the destination are protected (the C++ seed-only check could
    overwrite one). Local flood spreading is iterative, avoiding recursive stacks.
    """
    height, width = fixed.shape
    augmented = fixed.copy()
    known = np.isfinite(fixed)
    max_error = float(np.max(np.abs(surface[known] - fixed[known])))
    visited = [set(), set()]
    extrema = []
    for y in range(1, height - 1, step):
        for x in range(1, width - 1, step):
            if not support[y, x]:
                continue
            for direction in (-1, 1):
                seen = visited[direction == 1]
                cy, cx = y, x
                while 0 < cy < height - 1 and 0 < cx < width - 1:
                    if (cy, cx) in seen or not support[cy, cx]:
                        break
                    seen.add((cy, cx))
                    next_y, next_x = cy, cx
                    best = surface[cy, cx]
                    for dy, dx in _NEIGHBOURS:
                        candidate = surface[cy + dy, cx + dx]
                        if direction * (candidate - best) > 0.01:
                            best = candidate
                            next_y, next_x = cy + dy, cx + dx
                    if (next_y, next_x) == (cy, cx):
                        if not np.isfinite(augmented[cy, cx]):
                            # A perfectly flat plateau carries no extremum evidence.
                            neighbourhood = surface[cy - 1 : cy + 2, cx - 1 : cx + 2]
                            if np.ptp(neighbourhood) > 0.01:
                                augmented[cy, cx] = (
                                    surface[cy, cx] + direction * max_error * 2.5
                                )
                                extrema.append((cy, cx, direction))
                        break
                    cy, cx = next_y, next_x
    for y, x, _ in extrema:
        distances, differences = [], []
        for dy, dx in _NEIGHBOURS[:4]:
            for distance in range(1, max(height, width)):
                py, px = y + dy * distance, x + dx * distance
                if not (0 <= py < height and 0 <= px < width):
                    break
                if known[py, px]:
                    distances.append(distance)
                    differences.append(snapshot[py, px] - fixed[py, px])
                    break
        if not differences:
            continue
        correction = float(np.mean(differences))
        if abs(correction) <= 1e-8:
            continue
        radius = min(5, min(distances))
        stack, seen = [(y, x)], set()
        while stack:
            py, px = stack.pop()
            if (py, px) in seen or known[py, px]:
                continue
            seen.add((py, px))
            augmented[py, px] = snapshot[py, px] - correction
            for dy, dx in _NEIGHBOURS[:4]:
                ny, nx = py + dy, px + dx
                if (
                    0 <= ny < height
                    and 0 <= nx < width
                    and support[ny, nx]
                    and math.hypot(ny - y, nx - x) < radius - 1.1
                    and not np.isfinite(augmented[ny, nx])
                ):
                    stack.append((ny, nx))
    return augmented, {
        "minima": sum(d < 0 for _, _, d in extrema),
        "maxima": sum(d > 0 for _, _, d in extrema),
    }


def _bilinear(image, y, x):
    y = min(max(y, 0), image.shape[0] - 1)
    x = min(max(x, 0), image.shape[1] - 1)
    y0, x0 = int(y), int(x)
    y1, x1 = min(y0 + 1, image.shape[0] - 1), min(x0 + 1, image.shape[1] - 1)
    fy, fx = y - y0, x - x0
    return (1 - fy) * ((1 - fx) * image[y0, x0] + fx * image[y0, x1]) + fy * (
        (1 - fx) * image[y1, x0] + fx * image[y1, x1]
    )


def _gradient_constraints(surface, fixed, support, step=4):
    """Trace ascent/descent paths and interpolate between their height anchors."""
    known = np.isfinite(fixed)
    totals, counts = np.zeros_like(fixed), np.zeros(fixed.shape, dtype=np.int32)
    height, width = surface.shape
    accepted = 0
    for y in range(0, height, step):
        for x in range(0, width, step):
            if not known[y, x] or not support[y, x]:
                continue
            for direction in (-1, 1):
                px, py = float(x), float(y)
                ix, iy = x, y
                previous = surface[y, x]
                path, visited = [], {(iy, ix)}
                for distance in range(1, 1002):
                    current = _bilinear(surface, py, px)
                    if direction * (current - previous) < -1e-8:
                        break
                    previous = current
                    xp, xm = (
                        _bilinear(surface, py, px + 0.1),
                        _bilinear(surface, py, px - 0.1),
                    )
                    yp, ym = (
                        _bilinear(surface, py + 0.1, px),
                        _bilinear(surface, py - 0.1, px),
                    )
                    gx, gy = xp - xm, yp - ym
                    length = math.hypot(gx, gy)
                    if length <= 1e-6:
                        break
                    px += direction * gx / length
                    py += direction * gy / length
                    nx, ny = math.floor(px + 0.5), math.floor(py + 0.5)
                    if not (
                        0 < ny < height - 1 and 0 < nx < width - 1 and support[ny, nx]
                    ):
                        break
                    endpoints = [(ny, nx)]
                    if nx != ix and ny != iy:
                        endpoints += [(iy, nx), (ny, ix)]
                    endpoint = next(
                        ((ey, ex) for ey, ex in endpoints if known[ey, ex]), None
                    )
                    if endpoint is not None and (ny, nx) != (y, x):
                        if len(path) > 1:
                            end = fixed[endpoint]
                            for ry, rx, distance_from_start in path:
                                totals[ry, rx] += fixed[
                                    y, x
                                ] + distance_from_start / distance * (end - fixed[y, x])
                                counts[ry, rx] += 1
                            accepted += 1
                        break
                    if (ny, nx) != (iy, ix):
                        if (ny, nx) in visited:
                            break
                        visited.add((ny, nx))
                    iy, ix = ny, nx
                    if not known[iy, ix]:
                        path.append((iy, ix, distance))
    result = fixed.copy()
    traced = (counts > 0) & ~known
    result[traced] = totals[traced] / counts[traced]
    return result, {
        "accepted_paths": accepted,
        "added_constraint_pixels": int(traced.sum()),
    }


def interpolate(
    known_heights,
    unknown_ids=None,
    *,
    search_radius=120,
    refinement_iterations=100,
    gradient_step=4,
    extrema_step=40,
    max_error=1000.0,
    progress=None,
):
    """Run the full sequence; unknown heights are inferred once, then fixed.

    Unknown IDs may overlap known pixels to sample their heights, but cannot erase
    or change those known constraints. IDs need not be consecutive. No sea mask or
    implicit clipping of negative heights is applied.
    """
    fixed = np.asarray(known_heights, dtype=float).copy()
    if fixed.ndim != 2 or min(fixed.shape) < 3 or not np.isfinite(fixed).any():
        raise ValueError("Expected a 2-D grid with known heights and dimensions >= 3")
    if np.isinf(fixed).any():
        raise ValueError("Constraint heights must be finite or NaN")
    for value in (search_radius, refinement_iterations, gradient_step, extrema_step):
        if not isinstance(value, (int, np.integer)) or value <= 0:
            raise ValueError("Interpolation steps/radii must be positive integers")
    if not math.isfinite(max_error) or max_error <= 0:
        raise ValueError("Error scale must be positive and finite")
    unknown = (
        np.zeros(fixed.shape, dtype=np.int32)
        if unknown_ids is None
        else np.asarray(unknown_ids)
    )
    if (
        unknown.shape != fixed.shape
        or not np.issubdtype(unknown.dtype, np.integer)
        or (unknown < 0).any()
    ):
        raise ValueError(
            "Unknown contours require a matching nonnegative integer ID grid"
        )
    measured = np.isfinite(fixed)
    timings = {}
    started = time.monotonic()

    def stage(name):
        if progress:
            progress(name)

    stage("initial_directional_interpolation")
    initial = directional_fill(fixed, search_radius)
    support = np.isfinite(initial)
    working = _padding_surface(initial)
    timings["initial_seconds"] = time.monotonic() - started
    stage("unknown_height_assignment")
    estimates = []
    for contour_id in np.unique(unknown[unknown > 0]):
        mask = unknown == contour_id
        if not support[mask].all():
            raise ValueError(
                f"Unknown contour {contour_id} lacks initial-surface support"
            )
        values = initial[mask]
        height = float(values.mean())
        assign = mask & ~measured
        fixed[assign] = working[assign] = height
        estimates.append(
            {
                "unknown_id": int(contour_id),
                "estimated_elevation_m": height,
                "initial_min_m": float(values.min()),
                "initial_max_m": float(values.max()),
                "pixels": int(mask.sum()),
                "new_constraint_pixels": int(assign.sum()),
                "known_overlap_pixels": int((mask & measured).sum()),
            }
        )
    stage("initial_curvature_refinement")
    first, snapshot = _refine(working, fixed, refinement_iterations, max_error)
    stage("extrema_constraints")
    extrema, extrema_report = _extrema_constraints(
        first, fixed, snapshot, support, extrema_step
    )
    stage("gradient_constraints")
    gradients, gradient_report = _gradient_constraints(
        first, extrema, support, gradient_step
    )
    stage("final_directional_interpolation")
    second = _padding_surface(directional_fill(gradients, search_radius))
    stage("final_curvature_refinement")
    final, _ = _refine(second, gradients, refinement_iterations, max_error, final=True)
    assigned = np.isfinite(fixed)
    before_restore = float(np.max(np.abs(final[assigned] - fixed[assigned])))
    # Do not let smoothing replace a measured point or an established contour.
    final[assigned] = fixed[assigned]
    final[~support] = np.nan
    timings["total_seconds"] = time.monotonic() - started
    report = {
        "method": "imagico_python",
        "upstream_version": "0.3",
        "search_directions": 32,
        "search_radius_pixels": search_radius,
        "refinement_iterations": refinement_iterations,
        "gradient_step_pixels": gradient_step,
        "extrema_step_pixels": extrema_step,
        "known_constraint_pixels": int(measured.sum()),
        "assigned_unknown_contours": len(estimates),
        "unsupported_pixels": int((~support).sum()),
        "extrema": extrema_report,
        "gradients": gradient_report,
        "max_constraint_error_before_restoration_m": before_restore,
        "known_constraints_preserved": bool(
            np.array_equal(final[measured], np.asarray(known_heights)[measured])
        ),
        "unknown_height_assignment_passes": 1,
        "timing": timings,
    }
    return InterpolationResult(final, fixed, estimates, report)
