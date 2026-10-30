#!/usr/bin/env python3
"""Encode an already-resampled EGM2008 elevation tile as lossless Terrain-RGBA.

Usage: python terrain_rgb.py TILE_ELEVATIONS.tif OUTPUT.png|OUTPUT.webp
All pixels are opaque. Missing elevations display at zero metres; only the COGs
retain a distinguishable nodata value.
This step encodes one final-resolution raster. It does not reproject or build zooms.
"""
import argparse
import os
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

VALID_ALPHA = 255
# RGB for zero metres is a deterministic placeholder, not a missing height value.
NODATA_RGB = (1, 134, 160)


def encode_terrain_rgb(elevations, valid=None, *, round_bits=0):
    """Return opaque uint8 RGBA, replacing masked/nonfinite heights with zero.

    Source fill values must be masked BEFORE datum correction. Finite corrected
    heights equal to -255 or -1024 are eligible here, like any other valid height.
    round_bits erases low bits of the packed height code after numeric resampling,
    matching rio-rgbify's round-digits option. The nodata placeholder is never rounded.
    """
    if not isinstance(round_bits, (int, np.integer)) or not 0 <= round_bits <= 23:
        raise ValueError("round_bits must be an integer from 0 to 23")
    heights = np.ma.asarray(elevations, dtype=np.float64)
    if heights.ndim != 2:
        raise ValueError("Expected a two-dimensional elevation tile")
    data = np.ma.getdata(heights)
    usable = np.isfinite(data) & ~np.ma.getmaskarray(heights)
    if valid is not None:
        valid = np.asarray(valid, dtype=bool)
        if valid.shape != data.shape:
            raise ValueError("Validity mask must match the elevation tile shape")
        usable &= valid

    step = 2**round_bits
    q = np.rint(((data[usable] + 10000.0) / 0.1) / step) * step
    if np.any((q < 0) | (q > 16777215)):
        raise ValueError("Elevation is outside the Terrain-RGB encoding range")
    q = q.astype(np.uint32)
    rgba = np.empty((*data.shape, 4), dtype=np.uint8)
    rgba[..., :3] = NODATA_RGB
    rgba[..., 3] = VALID_ALPHA
    rgba[usable, 0] = (q >> 16) & 255
    rgba[usable, 1] = (q >> 8) & 255
    rgba[usable, 2] = q & 255
    rgba[usable, 3] = VALID_ALPHA
    return rgba


def write_terrain_rgb(elevations, output_path, valid=None, *, round_bits=0):
    """Write PNG or lossless WebP and verify that all four channels round-trip."""
    output_path = Path(output_path)
    formats = {".png": "PNG", ".webp": "WEBP"}
    image_format = formats.get(output_path.suffix.lower())
    if image_format is None:
        raise ValueError("Output must have a .png or .webp extension")
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite {output_path}")
    rgba = encode_terrain_rgb(elevations, valid, round_bits=round_bits)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=output_path.stem + ".", suffix=".part",
                                dir=output_path.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        options = {"lossless": True, "exact": True} if image_format == "WEBP" else {}
        Image.fromarray(rgba).save(temporary, format=image_format, **options)
        with Image.open(temporary) as check:
            if not np.array_equal(np.asarray(check.convert("RGBA")), rgba):
                raise ValueError("Encoded image did not preserve all four RGBA channels")
        if output_path.exists():
            raise FileExistsError(f"Output appeared during encoding: {output_path}")
        temporary.replace(output_path)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    import rasterio

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--round-bits", type=int, default=0,
                        help="Erase this many low height-code bits, after resampling (default: 0)")
    args = parser.parse_args()
    with rasterio.open(args.source) as src:
        if src.count != 1:
            raise ValueError("Expected a single-band, final-resolution elevation tile")
        heights = src.read(1, masked=True)
    write_terrain_rgb(heights, args.output, round_bits=args.round_bits)


if __name__ == "__main__":
    main()
