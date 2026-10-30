"""Small raster regressions for source fill, datum conversion, and RGBA storage."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import rasterio
from PIL import Image
from rasterio.transform import from_origin

import cartodem_to_egm2008 as converter
from terrain_rgb import encode_terrain_rgb, write_terrain_rgb


class NodataPipelineTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.geoid = self.root / "test_geoid.tif"
        with rasterio.open(self.geoid, "w", driver="GTiff", count=1,
                           width=8, height=8, dtype="float32", crs="EPSG:4979",
                           transform=from_origin(59, 36, 0.25, 0.25)) as dst:
            dst.write(np.full((8, 8), -50, dtype=np.float32), 1)
        self.pin = patch.object(converter, "GEOID_SHA256", converter.sha256(self.geoid))
        self.pin.start()
        self.addCleanup(self.pin.stop)

    def source(self, data, transform=None, nodata=None, mask=None):
        source = self.root / "source.tif"
        with rasterio.open(source, "w", driver="GTiff", count=1,
                           width=data.shape[1], height=data.shape[0], dtype="int16",
                           crs="EPSG:4326", nodata=nodata,
                           transform=transform or from_origin(60, 35, 1/3600, 1/3600)) as dst:
            dst.write(data.astype(np.int16), 1)
            if mask is not None:
                dst.write_mask(mask)
        return source

    def test_both_fill_values_survive_conversion_and_image_storage(self):
        # Cross a 512-pixel processing-window boundary as well as a height boundary.
        raw = np.full((2, 520), 100, dtype=np.int16)
        raw[0, :6] = [-255, -1024, 0, -66, -305, -1074]
        raw[0, 511:513] = [-255, -1024]
        mask = np.full(raw.shape, 255, dtype=np.uint8)
        mask[1, 0] = 0
        source = self.source(raw, mask=mask)
        output = self.root / "egm2008.tif"
        report = converter.convert(source, self.geoid, output)
        expected_invalid = np.isin(raw, [-255, -1024]) | (mask == 0)

        with rasterio.open(source) as original, rasterio.open(output) as corrected:
            heights = corrected.read(1, masked=True)
            np.testing.assert_array_equal(heights.mask, expected_invalid)
            self.assertTrue(np.isnan(corrected.nodata))
            self.assertEqual(corrected.transform, original.transform)
            self.assertEqual(corrected.crs, original.crs)
            np.testing.assert_allclose(heights.compressed(), (raw + 50)[~expected_invalid])
            self.assertEqual(report["source_nodata_pixels"], int(expected_invalid.sum()))
            self.assertEqual(report["source_nodata_values"], [-255, -1024])

        rgba = encode_terrain_rgb(heights)
        np.testing.assert_array_equal(rgba[..., 3], 255)
        np.testing.assert_array_equal(rgba[expected_invalid,:3], np.tile([1,134,160],(expected_invalid.sum(),1)))
        # Corrected heights that happen to equal source sentinels are still valid.
        self.assertEqual(heights[0, 4], -255)
        self.assertEqual(heights[0, 5], -1024)
        self.assertEqual(rgba[0, 4, 3], 255)
        self.assertEqual(rgba[0, 5, 3], 255)
        for suffix in ("png", "webp"):
            with self.subTest(format=suffix):
                path = self.root / f"tile.{suffix}"
                write_terrain_rgb(heights, path)
                with Image.open(path) as image:
                    np.testing.assert_array_equal(np.asarray(image.convert("RGBA")), rgba)

    def test_declared_nodata_is_also_preserved(self):
        source = self.source(np.array([[-9999, -255, -1024, 0]]), nodata=-9999)
        output = self.root / "egm2008.tif"
        converter.convert(source, self.geoid, output)
        with rasterio.open(output) as corrected:
            np.testing.assert_array_equal(corrected.read_masks(1), [[0, 0, 0, 255]])

    def test_missing_geoid_is_only_an_error_for_valid_source_pixels(self):
        outside = from_origin(80, 20, 1/3600, 1/3600)
        source = self.source(np.array([[-255, -1024]]), transform=outside)
        output = self.root / "empty.tif"
        report = converter.convert(source, self.geoid, output)
        self.assertEqual(report["valid_pixels"], 0)
        self.assertEqual(report["geoid_missing_pixels"], 0)
        with rasterio.open(output) as corrected:
            rgba = encode_terrain_rgb(corrected.read(1, masked=True))
            np.testing.assert_array_equal(rgba, [[[1,134,160,255],[1,134,160,255]]])

        source = self.source(np.array([[-1024, 0]]), transform=outside)
        failed_output = self.root / "missing_geoid.tif"
        with self.assertRaisesRegex(ValueError, "Geoid unavailable at 1"):
            converter.convert(source, self.geoid, failed_output)
        self.assertFalse(failed_output.exists())
        self.assertEqual(list(self.root.glob("*.part.tif")), [])

    def test_rgba_preserves_zero_negatives_and_explicit_validity(self):
        h = np.array([[0.0, -66.0, 12.34, np.nan, np.inf, 42.0]])
        valid = np.array([[True, True, True, True, True, False]])
        rgba = encode_terrain_rgb(h, valid)
        np.testing.assert_array_equal(rgba[..., 3], 255)
        rgb = rgba[..., :3].astype(np.int32)
        decoded = -10000 + .1 * (65536*rgb[..., 0] + 256*rgb[..., 1] + rgb[..., 2])
        np.testing.assert_allclose(decoded[0, :3], [0, -66, 12.3], atol=1e-6)
        np.testing.assert_array_equal(rgba[0, 3:], [[1, 134, 160, 255]] * 3)
        # Display tiles intentionally cannot distinguish measured zero from nodata.
        np.testing.assert_array_equal(rgba[0,0],rgba[0,3])
        with self.assertRaisesRegex(ValueError, "outside"):
            encode_terrain_rgb(np.array([[2_000_000.0]]))

    def water_mask(self, data, transform=None, nodata=0):
        path = self.root / "water.tif"
        with rasterio.open(path, "w", driver="GTiff", count=1,
                           width=data.shape[1], height=data.shape[0], dtype="uint8",
                           crs="EPSG:4326", nodata=nodata,
                           transform=transform or from_origin(60, 35, 1/3600, 1/3600)) as dst:
            dst.write(data.astype(np.uint8), 1)
        return path

    def test_coastal_mask_keeps_land_lakes_and_rivers_with_finite_nodata(self):
        raw = np.array([[0, 10, 20, 30, -255, -1024, -60]])
        source = self.source(raw)
        # Product nodata metadata can say zero, although class zero is valid land.
        water = self.water_mask(np.array([[0, 1, 2, 3, 0, 1, 0]]))
        output = self.root / "coastal.tif"
        report = converter.convert(source, self.geoid, output, water_mask_path=water,
                                   output_nodata=-255)
        with rasterio.open(output) as src:
            self.assertEqual(src.nodata, -255)
            np.testing.assert_array_equal(src.read(1), [[50, -255, 70, 80, -255, -255, -10]])
            h = src.read(1, masked=True)
            np.testing.assert_array_equal(h.mask, [[False, True, False, False, True, True, False]])
        self.assertEqual(report['source_fill_counts'], {'-255': 1, '-1024': 1})
        self.assertEqual(report['newly_masked_ocean_pixels'], 1)
        self.assertEqual(report['source_nodata_on_glo30_land_pixels'], 1)
        self.assertEqual(report['valid_pixels'], 4)
        rgba = encode_terrain_rgb(h, round_bits=11)
        np.testing.assert_array_equal(rgba[h.mask], [[1,134,160,255]] * 3)

    def test_finite_nodata_collision_fails_without_publishing_sheet(self):
        source = self.source(np.array([[-305, -255]]))
        output = self.root / 'collision.tif'
        with self.assertRaisesRegex(ValueError, 'collides'):
            converter.convert(source, self.geoid, output, output_nodata=-255)
        self.assertFalse(output.exists())

    def test_water_mask_missing_coverage_fails(self):
        source = self.source(np.array([[10, 20, 30]]))
        water = self.water_mask(np.array([[0, 1]]))
        output = self.root / 'incomplete_mask.tif'
        with self.assertRaisesRegex(ValueError, 'missing coverage'):
            converter.convert(source, self.geoid, output, water_mask_path=water,
                              output_nodata=-255)
        self.assertFalse(output.exists())

    def test_water_mask_alignment_uses_pixel_centres(self):
        source = self.source(np.array([[10, 20, 30, 40]]))
        # Half as many mask pixels, at twice the source spacing.
        water = self.water_mask(np.array([[1, 0]]),
                                transform=from_origin(60,35,2/3600,2/3600))
        output = self.root / 'aligned_mask.tif'
        converter.convert(source, self.geoid, output, water_mask_path=water, output_nodata=-255)
        with rasterio.open(output) as src:
            np.testing.assert_array_equal(src.read(1), [[-255,-255,80,90]])

    def test_zoom_precision_rounds_heights_and_preserves_nodata_bytes(self):
        heights = np.array([[0.0, 12.23, -66.0, 1000.0, np.nan]])
        for bits in range(4, 12):
            with self.subTest(round_bits=bits):
                rgba = encode_terrain_rgb(heights, round_bits=bits)
                rgb = rgba[0, :4, :3].astype(np.int32)
                codes = rgb[:, 0]*65536 + rgb[:, 1]*256 + rgb[:, 2]
                np.testing.assert_array_equal(codes % 2**bits, 0)
                errors = np.abs(-10000 + .1*codes - heights[0, :4])
                self.assertLessEqual(float(errors.max()), .05*2**bits + 1e-9)
                np.testing.assert_array_equal(rgba[0, 4], [1, 134, 160, 255])
                np.testing.assert_array_equal(rgba[0, :4, 3], 255)
                if bits >= 8:
                    np.testing.assert_array_equal(rgba[0, :4, 2], 0)
                output = self.root / f"rounded_{bits}.webp"
                write_terrain_rgb(heights, output, round_bits=bits)
                with Image.open(output) as image:
                    np.testing.assert_array_equal(np.asarray(image.convert("RGBA")), rgba)
        # Quantization uses the offset packed code, not a zero-anchored metre grid.
        np.testing.assert_array_equal(encode_terrain_rgb(heights, round_bits=8)[0, 0],
                                      [1, 135, 0, 255])
        for invalid_bits in (-1, 24, 1.5):
            with self.assertRaisesRegex(ValueError, "round_bits"):
                encode_terrain_rgb(heights, round_bits=invalid_bits)


if __name__ == "__main__":
    unittest.main()
