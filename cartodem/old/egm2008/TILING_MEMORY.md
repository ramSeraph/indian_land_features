# TerrainRGB memory correction, 25 October 2026

All 52,807 river components, 543 candidate COGs and the final mosaic completed.
The first TerrainRGB run was killed for insufficient memory during zoom 5,
after committing four tiles. The kernel recorded 1,649,996 KiB of anonymous
resident memory for the killed process. Completed COGs are reused unchanged.

## Cause and correction

The 136,800 by 111,600 native-pixel mosaic has no mosaic-level overviews. A
low-zoom warped VRT block can require a very large native source window.
GDAL documents that `WarpRegionToBuffer` processes a region in one chunk,
whereas `ChunkAndWarpImage` subdivides it according to the warp memory budget:
[GDAL warp operation reference](https://gdal.org/en/stable/doxygen/classGDALWarpOperation.html).
A 32 MiB warp setting alone did not bound the previous warped VRT reads.

`production_tiles.Sampler` now uses `rasterio.warp.reproject` with a 32 MiB
warp budget, a 64 MiB GDAL block cache and one thread. GDAL subdivides the
source reads while producing each 512-pixel tile. Every tile transform is
derived from the same zoom-level grid. The existing transformation tolerance
of 0.000001 source pixels is retained.

Numeric resampling remains average at zooms 5–11 and bilinear at zoom 12.
The separate native-resolution validity mosaic still uses nearest pixel-centre
sampling. Missing pixels keep `(1, 134, 160, 254)` and valid pixels alpha 255;
numeric resampling precedes `round_bits = 16-z` quantization and lossless WebP.
No COG, river fit, mask or nodata policy is changed.

## Validation

- The bounded reproduction of the original reader exceeded the 1,200 MiB RSS
  diagnostic threshold and was terminated after 32.08 seconds; sampled peak
  RSS was 1,211.53 MiB. This sampler is a diagnostic guard, not a hard limit.
- The initial chunked-reader trial processed all 16 zoom-5 tiles in 252.76
  seconds, with 438.65 MiB peak observed RSS. Its four overlapping tiles are
  byte-identical to the previously committed WebP tiles, including RGBA and
  validity counts. This trial used Rasterio's default transform tolerance;
  the final implementation restores the preceding 0.000001-pixel tolerance.
- All 63 tests passed with the final implementation. A new test spans both
  source-sheet and XYZ tile boundaries, a varying elevation surface and a
  nodata hole. At zooms 9 and 12, validity is identical to the preceding
  whole-grid warped VRT reference and valid heights agree within 0.0001 m.
  Existing tests cover resuming tile batches and exact RGBA in PMTiles.

Portable measurements and the restart record are under
[diagnostics/tile-memory-v1](diagnostics/tile-memory-v1/). To repeat the full
zoom-5 read/encode trial against an existing completed mosaic:

```bash
python tile_memory_trial.py --workdir /path/to/build --label chunked-check
python -m unittest discover -v
```

Use a new trial label to preserve earlier measurements. The trial writes logs
and hashes under `reports/tile-memory-<label>/`; it does not alter production
tile databases.

## Restart

Preserve the failed status and the four-tile `tiles/final/` checkpoint under
the archive paths recorded in `diagnostics/tile-memory-v1/restart.json`.
Start a fresh tile database because the implementation and recipe have changed;
do not rewrite old recipe hashes. Inputs, COGs, masks, mosaic and river
checkpoints remain in place. Resume with:

```bash
python run_build.py --workdir /path/to/build --river-source glo30 \
  --resume-from tiles --detach
```

The build remains a candidate pending the existing river-gradient review.
The runner creates output files; it does not publish a release.

## Extending the completed pyramid to zoom 0

The later z0–4 extension exposed another allocation failure on the world tile:
GDAL requested 2,175,120,004 bytes during the numeric warp, despite the 32 MiB
warp budget. The attempt failed before committing a tile. The existing z5–12
archives remained intact. The first attempt is recorded in
`reports/lowzoom-final-first-attempt.json` and `logs/lowzoom-final.log`.

For z0–4, the sampler first reads the native validity mask as before. It then
warps heights only into the bounding rectangle of valid destination pixels.
This removes the world tile's empty margins from the numeric warp while
retaining the same destination pixel edges, native source data and numeric
average. The rest of the output remains nodata. Zooms 5–12 keep their previous
sampling path. This avoids introducing a separate approximate overview
surface or resizing packed RGB values.

The focused production, pipeline and launcher checks passed (18 tests). The
new world-tile test compares valid heights against an uncropped numeric warp
within 0.001 m. Extension tests also check archive preservation, unchanged
native nodata treatment, quantization, repeated invocation and rejection of
changed source rasters. Zooms 0–4 retain the z5 quantization step of 204.8 m:
`round_bits = 16 - max(z, 5)`.

The extension command and preserved baseline are documented in [README.md](README.md).
`reports/lowzoom-final.json` records completion and peak resident memory;
`reports/source-lowzoom-cropped.tar.gz` preserves the source used for this run.

The retry completed all 12 added tiles and their archive on 25 October 2026.
Each zoom took 238–269 seconds; peak resident memory was 854,796 KiB
(834.76 MiB). Including input/archive verification, the extension took
26 minutes 15 seconds. The five existing archive records and checksums were
preserved. Completion records are copied to [records/final](records/final/).
