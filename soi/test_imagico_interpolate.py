# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "scipy>=1.14", "pytest>=8"]
# ///
import numpy as np
import pytest

from imagico_interpolate import directional_fill, interpolate


def plane():
    heights = np.full((41, 41), np.nan)
    heights[::20, :] = np.arange(0, 41, 20)[:, None]
    return heights


def test_directional_search_preserves_zero_and_negative_constraints():
    known = plane() - 20
    output = directional_fill(known)
    assert np.array_equal(output[np.isfinite(known)], known[np.isfinite(known)])
    assert np.allclose(output[10, 8:-8], -10)
    assert np.allclose(output[30, 8:-8], 10)


def test_unknown_contour_is_assigned_once_and_used_in_final_surface():
    known = plane()
    snapshot = known.copy()
    unknown = np.zeros(known.shape, dtype=np.int32)
    unknown[10, 8:-8] = 17
    unknown[30, 8:-8] = 200
    result = interpolate(known, unknown)
    assert np.array_equal(known, snapshot, equal_nan=True)
    assert [r["unknown_id"] for r in result.unknown_heights] == [17, 200]
    assert [
        r["estimated_elevation_m"] for r in result.unknown_heights
    ] == pytest.approx([10, 30])
    assert np.all(result.elevation[10, 8:-8] == 10)
    assert np.all(result.elevation[30, 8:-8] == 30)
    assert result.report["unknown_height_assignment_passes"] == 1
    assert result.report["known_constraints_preserved"]
    assert result.report["gradients"]["added_constraint_pixels"] > 0
    assert np.max(np.abs(result.elevation[:, 10:-10] - np.arange(41)[:, None])) < 2


def test_unknown_height_does_not_replace_a_known_point_at_an_overlap():
    known = plane()
    known[10, 20] = 11.5
    unknown = np.zeros(known.shape, dtype=np.int32)
    unknown[10, 8:-8] = 1
    result = interpolate(known, unknown)
    assert result.elevation[10, 20] == 11.5
    assert result.constraints[10, 20] == 11.5
    assert result.unknown_heights[0]["known_overlap_pixels"] == 1
    assert result.report["known_constraints_preserved"]


def test_surface_is_invariant_to_a_vertical_offset_including_below_sea_level():
    known = plane()
    unknown = np.zeros(known.shape, dtype=np.int32)
    unknown[10, 8:-8] = 1
    a = interpolate(known, unknown)
    b = interpolate(known - 25, unknown)
    assert np.allclose(a.elevation - 25, b.elevation, atol=1e-7)
    assert b.unknown_heights[0]["estimated_elevation_m"] == pytest.approx(-15)


def test_nested_contours_form_a_summit_and_a_symmetric_depression():
    known = np.full((81, 81), np.nan)
    y, x = np.mgrid[:81, :81]
    distance = np.hypot(x - 40, y - 40)
    for radius, level in [(30, 100), (20, 120), (10, 140)]:
        known[abs(distance - radius) < 0.55] = level
    hill = interpolate(known, extrema_step=10)
    depression = interpolate(-known, extrema_step=10)
    assert hill.report["extrema"]["maxima"] > 0
    assert depression.report["extrema"]["minima"] > 0
    assert hill.elevation[40, 40] > 140
    assert np.allclose(hill.elevation, -depression.elevation, equal_nan=True)
    assert hill.report["known_constraints_preserved"]
    assert depression.report["known_constraints_preserved"]


def test_sparse_support_is_not_filled_with_fake_zero_heights():
    known = np.full((31, 31), np.nan)
    known[15, 15] = 100
    result = directional_fill(known, radius=3)
    assert np.isnan(result[0, 0])
    assert result[15, 15] == 100
    unknown = np.zeros(known.shape, dtype=np.int32)
    unknown[1, 1] = 1
    with pytest.raises(ValueError, match="lacks initial-surface support"):
        interpolate(known, unknown, search_radius=3)


@pytest.mark.parametrize("bad", [-1, 0, 2.5])
def test_invalid_parameters_are_rejected(bad):
    with pytest.raises(ValueError):
        interpolate(plane(), refinement_iterations=bad)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
