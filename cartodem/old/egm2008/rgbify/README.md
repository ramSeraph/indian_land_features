# TerrainRGB using the original command-line workflow

This is an evaluated alternative, retained with its pilot results. The selected
workflow keeps the existing numeric tiler and exports its encoded databases to
MBTiles before running the pmtiles-mosaic CLI; see `../TILING_SETUP.md`.

This workflow starts with the completed EGM2008 COGs, including ocean clipping,
river treatment, and NoData=-255. It follows `../../steps.sh`: `gdalbuildvrt`,
`gdalwarp` to Web Mercator with cubic-spline resampling and destination nodata 0,
then `rio rgbify` producing WebP MBTiles for each zoom.

The zoom range is extended to 0–12. The original quantization is retained:
`round_bits=16-max(z,5)`. `pmtiles-mosaic` accepts the thirteen MBTiles files
as separate `--from-source` arguments and creates geographic PMTiles partitions.
There is no custom PMTiles writer. Images are 512 pixels and lossless WebP.

The current script uses stock rio-rgbify. Its quantization also rounds the zero
nodata placeholder: at z0–6 it becomes 35.2 m, while at z11–12 it remains 0 m.
The full dataset was not rebuilt with this alternative, and no rio-rgbify patch
was applied.

## Install

Install `gdal-bin`, `git`, and `uv`. Create separate environments because
`pmtiles-mosaic` pins its PMTiles library dependency independently:

```bash
sudo apt-get update
sudo apt-get install -y gdal-bin git
uv venv --python 3.14.3 .venv-rgbify
uv pip sync --python .venv-rgbify/bin/python requirements-rgbify.txt
uv venv --python 3.14.3 .venv-mosaic
uv pip sync --python .venv-mosaic/bin/python requirements-mosaic.txt
```

The Pillow wheel includes WebP support. Verify it before running:

```bash
.venv-rgbify/bin/python -c 'from PIL import features; assert features.check("webp")'
```

## Run

```bash
RIO_BIN="$PWD/.venv-rgbify/bin/rio" \
PARTITION_BIN="$PWD/.venv-mosaic/bin/partition" \
bash run.sh /path/to/final/cogs /path/to/output 4
```

COGs are read-only. Completed MBTiles files are checked and reused on restart.
Incomplete MBTiles files cause a stop. The coverage check is necessary because
the upstream tiler logs individual worker failures without always failing the
whole command. Keep the MBTiles until all PMTiles payloads and mosaic coverage
have been verified. Uploading is a separate step and does not change release notes.
