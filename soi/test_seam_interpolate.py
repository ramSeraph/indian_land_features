# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "scipy>=1.14", "pytest>=8"]
# ///
import numpy as np
import pytest

from seam_interpolate import reconcile_patch


def test_local_patch_preserves_source_heights_and_unchanged_surroundings():
    known = np.full((41, 51), np.nan)
    known[:, 5] = 50
    known[:, 45] = 50
    known[20, 25] = 50
    prior = np.full(known.shape, 48.0)
    prior[:, 25:] = 52
    prior[np.isfinite(known)] = known[np.isfinite(known)]
    core = np.zeros(known.shape, bool)
    core[2:-2, 15:36] = True
    patched, report = reconcile_patch(
        known, np.zeros(known.shape, dtype=int), prior, core
    )
    np.testing.assert_array_equal(patched[~core], prior[~core])
    assert patched[20, 25] == 50
    assert np.isfinite(patched[core]).all()
    assert np.any(patched[core] != prior[core])
    assert report["measured_constraints_preserved"]
    assert report["outside_core_unchanged"]
    assert report["artificial_boundary_cells"] > 0


def test_unknown_group_crossing_replacement_boundary_requires_expansion():
    known = np.full((21, 21), np.nan)
    unknown = np.zeros(known.shape, dtype=int)
    unknown[10, 3:18] = 1
    core = np.zeros(known.shape, bool)
    core[5:16, 5:16] = True
    with pytest.raises(ValueError, match="enlarge"):
        reconcile_patch(known, unknown, np.full(known.shape, 50.0), core)


def test_complete_unknown_groups_are_reestimated_and_outside_groups_retained():
    known = np.full((21, 21), np.nan)
    known[5, 5] = 50
    unknown = np.zeros(known.shape, dtype=int)
    unknown[10, 9:12] = 1
    unknown[1, 1:4] = 2
    core = np.zeros(known.shape, bool)
    core[5:16, 5:16] = True
    prior = np.full(known.shape, 50.0)
    patched, report = reconcile_patch(known, unknown, prior, core)
    assert report["inherited_unknown_context_cells"] == 3
    estimates = report["unknown_height_estimates"]
    assert len(estimates) == 1 and estimates[0]["unknown_id"] == 1
    assert estimates[0]["estimated_elevation_m"] == pytest.approx(50)
    np.testing.assert_allclose(patched, 50)


def test_four_quadrant_corner_is_solved_together():
    # A synthetic four-sheet join with distinct independent-sheet offsets.
    rows, cols = np.indices((61, 61))
    prior = 100.0 + np.where(cols < 30, -3, 3) + np.where(rows < 30, -2, 2)
    known = np.full(prior.shape, np.nan)
    known[10, 10] = known[10, 50] = known[50, 10] = known[50, 50] = 100
    prior[np.isfinite(known)] = known[np.isfinite(known)]
    core = (abs(rows - 30) < 10) | (abs(cols - 30) < 10)
    core[:2, :] = core[-2:, :] = False
    core[:, :2] = core[:, -2:] = False
    patched, report = reconcile_patch(
        known, np.zeros(prior.shape, dtype=int), prior, core
    )
    before_steps = np.concatenate(
        [prior[29, 20:41] - prior[30, 20:41], prior[20:41, 29] - prior[20:41, 30]]
    )
    after_steps = np.concatenate(
        [
            patched[29, 20:41] - patched[30, 20:41],
            patched[20:41, 29] - patched[20:41, 30],
        ]
    )
    assert np.mean(abs(after_steps)) < np.mean(abs(before_steps)) / 5
    np.testing.assert_array_equal(patched[~core], prior[~core])
    assert report["outside_core_unchanged"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
