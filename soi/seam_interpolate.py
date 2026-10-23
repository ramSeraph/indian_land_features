# SPDX-License-Identifier: GPL-3.0-or-later
"""Local Imagico reconciliation with explicit fixed outer boundary values.

Inputs are constraints on one common metric grid, never a blend of final DEMs.
The prior surface supplies artificial outer-boundary conditions only. Keep
original source constraints and these boundary conditions distinct in provenance.
"""

import numpy as np
from scipy.ndimage import binary_erosion

from imagico_interpolate import interpolate


def reconcile_patch(known, unknown, prior, core, *, frame_pixels=2, progress=None):
    known = np.asarray(known, dtype=float)
    prior = np.asarray(prior, dtype=float)
    unknown = np.asarray(unknown)
    core = np.asarray(core, dtype=bool)
    if not (known.shape == prior.shape == unknown.shape == core.shape):
        raise ValueError("Patch inputs must share a grid")
    if frame_pixels < 1 or not core.any():
        raise ValueError(
            "Expected a nonempty replacement core and positive frame width"
        )
    for group in np.unique(unknown[core & (unknown > 0)]):
        if np.any((unknown == group) & ~core):
            raise ValueError(
                f"Unknown-height group {group} crosses the replacement boundary; enlarge the region"
            )
    frame = core & ~binary_erosion(core, iterations=frame_pixels, border_value=0)
    boundary = frame & np.isfinite(prior) & ~np.isfinite(known) & (unknown == 0)
    fixed = known.copy()
    # Re-estimate complete unknown groups inside the replacement region. Outside
    # groups retain their first-pass estimates as labelled context, rather than
    # computing new heights that would not be written back to those contours.
    inherited = ~core & (unknown > 0) & np.isfinite(prior) & ~np.isfinite(known)
    fixed[inherited] = prior[inherited]
    unknown = unknown.copy()
    unknown[~core] = 0
    fixed[boundary] = prior[boundary]
    result = interpolate(fixed, unknown, progress=progress)
    if np.any(core & np.isfinite(prior) & ~np.isfinite(result.elevation)):
        raise ValueError("Local solve loses supported cells; enlarge the context")
    patched = prior.copy()
    patched[core] = result.elevation[core]
    measured = core & np.isfinite(known)
    if not np.array_equal(patched[measured], known[measured]):
        raise ValueError("Local solve changed a measured height constraint")
    if not np.array_equal(patched[~core], prior[~core], equal_nan=True):
        raise ValueError("Local solve changed cells outside the replacement region")
    report = {
        "method": "imagico_local_constraints_with_fixed_prior_boundary",
        "replacement_cells": int(core.sum()),
        "artificial_boundary_cells": int(boundary.sum()),
        "inherited_unknown_context_cells": int(inherited.sum()),
        "boundary_frame_pixels": frame_pixels,
        "measured_constraint_pixels_in_core": int(measured.sum()),
        "measured_constraints_preserved": True,
        "outside_core_unchanged": True,
        "unknown_height_estimates": result.unknown_heights,
        "engine_report": result.report,
        "limitations": [
            "Fixed outer-boundary values are inherited model estimates, not new measurements.",
            "This does not repair inconsistent source heights or establish equality to country-wide interpolation.",
        ],
    }
    return patched, report
