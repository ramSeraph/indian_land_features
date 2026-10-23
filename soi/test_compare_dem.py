# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "rasterio>=1.4", "scipy>=1.14", "pyproj>=3.6", "matplotlib>=3.9", "pytest>=8"]
# ///
import numpy as np
import pytest
from rasterio.transform import Affine

from compare_dem import align_bounds, sample_bilinear, stats


def test_window_preserves_native_pixel_origin_and_encloses_requested_area():
    transform = Affine(1 / 3600, 0, 84 - 1 / 7200, 0, -1 / 3600, 21 + 1 / 7200)
    requested = (
        84.74485533024225,
        20.245515221656664,
        85.00490498990867,
        20.504641972054507,
    )
    left, bottom, right, top = align_bounds(requested, transform)
    assert left <= requested[0] and bottom <= requested[1]
    assert right >= requested[2] and top >= requested[3]
    col0, row0 = ~transform * (left, top)
    col1, row1 = ~transform * (right, bottom)
    np.testing.assert_allclose(
        [col0, row0, col1, row1], np.rint([col0, row0, col1, row1]), atol=1e-8
    )


def test_bilinear_sampling_at_pixel_centres_and_between_pixels():
    transform = Affine(30, 0, 500000, 0, -30, 2200000)
    rows, cols = np.indices((4, 5))
    x, y = transform * (cols + 0.5, rows + 0.5)
    values = 0.2 * (x - 500000) - 0.1 * (y - 2200000) + 100
    points_col = np.array([0.5, 2.25, 4.5])
    points_row = np.array([0.5, 1.75, 3.5])
    sx, sy = transform * (points_col, points_row)
    actual = sample_bilinear(values, transform, sx, sy)
    expected = 0.2 * (sx - 500000) - 0.1 * (sy - 2200000) + 100
    np.testing.assert_allclose(actual, expected, atol=1e-8)
    assert np.isnan(
        sample_bilinear(values, transform, np.array([499000]), np.array([2200000]))[0]
    )


def test_difference_sign_statistics_and_invalid_samples():
    result = stats([2, -2, np.nan, 0, np.inf])
    assert result["count"] == 3
    assert result["mean_m"] == 0
    assert result["mae_m"] == pytest.approx(4 / 3)
    assert result["rmse_m"] == pytest.approx(np.sqrt(8 / 3))


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
