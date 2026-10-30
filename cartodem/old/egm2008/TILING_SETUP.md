# Rebuild display tiles from existing COGs on ARM Linux

The current publication path retains the numeric tiler, exports its per-zoom
SQLite databases to standard MBTiles, and runs the existing `pmtiles-mosaic`
partition CLI with thirteen sources. Use `--tiles-only` when running
`rebuild_tiles.py`; the earlier custom PMTiles writer is not used for publication.
See the commands in "Standard MBTiles and PMTiles packaging" below.

This workflow reuses the final EGM2008 COGs after ocean clipping and river
treatment. It does not download source DEMs or repeat datum conversion, routing
or river fitting. COGs retain native Float32 values and declared NoData=-255.
TerrainRGB is a display product: missing pixels become exactly 0 m with full
opacity. There is no alpha-254 or 0.1 m sentinel. Keep the existing quantization,
including 1.6 m increments at z12 (maximum rounding error 0.8 m).

## Software

Use 64-bit ARM Linux with glibc >= 2.28, preferably Ubuntu 22.04/24.04 or Debian
12. Install the small bootstrap prerequisites if the image lacks them:

```bash
sudo apt-get update
sudo apt-get install -y ca-certificates curl tar gzip coreutils
bash egm2008/install_tiles_arm64.sh /path/to/tiling-tools
```

The installer checks the architecture and downloads checksum-pinned Linux ARM64
binaries for uv 0.13.6 and GitHub CLI 2.103.2. uv installs Python 3.14.3 and a
fresh virtual environment with [requirements-tiles.txt](requirements-tiles.txt):
NumPy, Rasterio, Pillow, the Python PMTiles library and their pinned dependencies.
Binary wheels are required, so it never silently compiles GDAL from source.
Rasterio includes GDAL/PROJ and the GTiff/COG/VRT drivers; Pillow includes WebP.
The installer also creates `mosaic-venv` with `pmtiles-mosaic==0.2.0` and its
dependencies from [requirements-mosaic.txt](requirements-mosaic.txt).
GRASS, WhiteboxTools, Rust, Node, rio-rgbify, the PMTiles Go CLI and a separate
GDAL command-line installation are not needed for this workflow.

Everything is installed under the selected tools directory, including the
Python interpreter and download cache. No shell profile or system Python is
changed. The installer verifies lossless WebP, opaque-zero pixels, GDAL drivers
and PMTiles indexing; `software.json` records installed versions. Authentication
is separate and is not stored in scripts:

```bash
/path/to/tiling-tools/bin/gh auth login --hostname github.com --git-protocol https --web
/path/to/tiling-tools/bin/gh auth status
```

Use a fresh environment rather than moving the old virtual environment between
OS installations. See [uv installation](https://docs.astral.sh/uv/getting-started/installation/)
and [Rasterio installation](https://rasterio.readthedocs.io/en/stable/installation.html).

## CPU, memory, storage and time

| Resource | Guidance for the tiling-only rebuild |
| --- | --- |
| CPU | Choose 4 ARM64 vCPUs and `--workers 4` for a parallel rebuild. One worker still works on one CPU. |
| RAM | 8 GB is a reasonable starting allocation for four workers; 16 GB gives substantial headroom. For one worker, 2 GB worked previously and 4 GB is comfortable. |
| More CPU/RAM | `--workers N` runs N isolated tile processes at z5–12. Allow roughly 2 GB of RAM per worker as a conservative initial allocation and monitor total memory. This is a planning estimate, not a measured four-worker peak. |
| GPU | Not used. |
| Free disk | Allow at least 30 GiB beyond the existing COGs and anything retained from the old build. |
| Disk performance | Prefer SSD/block storage with good read throughput; the COG mosaic is repeatedly sampled. |
| Time | The single-worker baseline was about 8 hours. Parallel throughput depends on CPU and disk bandwidth; no linear speedup is promised. |

Recorded measurements from the previous build:

- All zooms 0–12: 26,939.61 seconds (7 h 29 min) of tile generation, for
  227,426 tiles. z12 accounted for 18,206.78 seconds (5 h 3 min).
- The z0–4 extension peaked at 854,796 KiB RSS (834.76 MiB). An earlier z5
  reader/encoder trial peaked at 438.65 MiB. These are measured stages, not a
  hard memory limit for every future dataset or GDAL version.
- The completed PMTiles occupied about 7.27 GB. Retaining the encoded SQLite
  databases and their standard MBTiles exports requires two further copies
  of roughly that size. The partitioner also uses temporary disk space.
  Allow 30 GiB for these copies, temporary files, masks and software.

See [the timing record](records/final/final-tiles.json),
[the low-zoom memory record](records/final/lowzoom-final.json) and
[the memory implementation notes](TILING_MEMORY.md). The old numbers are estimates
for later rebuilds. The opaque four-worker build completed in about 2 h 15 min,
including its original packaging stage.

## Rebuild

Keep both the COG directory and its original `rivers/final/cogs.json` manifest
from the completed build. Use a separate output workspace so new tiles never
mix with the previous alpha-254 archives:

```bash
/path/to/tiling-tools/venv/bin/python egm2008/rebuild_tiles.py \
  --workdir /path/to/new-build \
  --cog-dir /path/to/old-build/output/final/cogs \
  --cog-manifest /path/to/old-build/rivers/final/cogs.json \
  --workers 4 \
  --tiles-only \
  --detach
```

The default expects all 543 sheets. `--expected-sheets N --min-zoom Z
--max-zoom Z` supports a deliberately smaller pilot with its own manifest.

The runner verifies every source file against its original SHA-256 and byte
size, creates a VRT and native validity masks, generates opaque 512-pixel
lossless WebP tiles at z0–12, and leaves them in per-zoom SQLite databases.
The separate packaging commands below create MBTiles and PMTiles.
Numeric resampling and quantization are unchanged. COGs are opened
read-only. Genuine zero heights and missing heights are intentionally
indistinguishable in the display tiles. Consumers use ordinary `raster-dem`
decoding with `encoding: "mapbox"` and `tileSize: 512`.

Progress and output paths, relative to the new workspace:

- `logs/rebuild-display.log`: stage and tile progress.
- `reports/rebuild-display.json`: process status, final totals, parent peak RSS
  and child `ru_maxrss` (the largest child-process peak on Linux, not aggregate RAM).
- `reports/display-tile-status.json`: current zoom and completed tile count.
- `tiles/display/`: resumable per-zoom databases.
- `output/display/`: tile metadata, COG manifest, checksums, runtime versions,
  rebuild notes and source archive. With `--tiles-only`, PMTiles are created
  separately in the directory passed to `partition_mbtiles.sh`.

Rerun the same command after an interruption. Committed tile columns are reused;
changed input checksums, code/encoding recipes or zoom ranges are rejected.
The worker count can change on restart. Tasks contain at most 32 tile rows,
with at most twice the worker count submitted at once. Each process owns its
GDAL datasets and uses one GDAL thread; only the parent writes SQLite. Partial
columns from an interruption are regenerated. Zooms 0–4 and final packaging
remain serial. Parallelism is independent of the alpha/nodata encoding.
No uploads are performed by the installer or the rebuild command.

## Standard MBTiles and PMTiles packaging

The per-zoom databases contain encoded WebP bytes. Export changes their schema
and converts XYZ y coordinates to MBTiles TMS row numbers. It does not change
image bytes or rerun DEM processing. Every exported tile is checked against its
original SHA-256; all thirteen resulting MBTiles sources are passed directly to
the existing partition CLI:

```bash
python egm2008/export_mbtiles.py --workdir /path/to/build \
  --outdir /path/to/build/publication/mbtiles

bash egm2008/partition_mbtiles.sh /path/to/build/publication/mbtiles \
  /path/to/build/publication/pmtiles /path/to/mosaic-env/bin/partition

/path/to/mosaic-env/bin/python egm2008/verify_partitioned.py \
  --workdir /path/to/build --pmtiles-dir /path/to/build/publication/pmtiles
```

Use the installer's `mosaic-venv` for `partition` and `verify_partitioned.py`.
The verifier records its complete dependency versions in `requirements-mosaic.txt`. The standard
partitioner receives a 2,000,000,000-byte file limit and chooses geographic
partitions. The verifier checks every packaged tile against the encoded-tile
databases. Validate the generated mosaic with `verify_mosaic.mjs` and the same
PMTiles mosaic reader used by the serving application before uploading.

For replacement publishing, pass this directory as `--pmtiles-dir` to
`replace_release.py` during preparation. It retains the mosaic generated by the
standard tool and maps the archive names to the existing release prefix.

The previous preview uses the historical alpha-254 archives. A height-difference
view of the new display tiles must obtain coverage from the COGs or their
validity raster; it must not infer missing data from RGB zero.

## Publication

The selected workflow replaces the existing release as described below.
Individual `.tif` assets keep the COGs directly usable; the PMTiles archives
are already partitioned.
GitHub permits up to 1,000 assets per release, each under 2 GiB, so 543 COGs
plus manifests fit within one release. The upload helper checks those limits.
See [GitHub release limits](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases#storage-and-bandwidth-quotas)
and [gh release upload](https://cli.github.com/manual/gh_release_upload).

## Replace assets in an existing release

`replace_release.py` implements an explicitly requested replacement while
preserving the release title, notes, tag, and publication flags. It first checks
all COG and PMTiles sizes and SHA-256 digests, snapshots the release inventory,
and prepares compatible manifests. COG names are preserved; PMTiles are mapped
to the supplied prefix with `-part0000.pmtiles`, `-part0001.pmtiles`, and so on.
The corresponding `.mosaic.json`, `.partition_info.json`, and `listing_files.csv`
are regenerated from the new files. The original build files are not renamed.

```bash
python egm2008/replace_release.py \
  --workdir /path/to/new-build \
  --pmtiles-dir /path/to/new-build/publication/pmtiles \
  --repo OWNER/REPOSITORY --tag RELEASE_TAG \
  --prefix EXISTING_TERRAINRGB_PREFIX \
  --gh /path/to/tiling-tools/bin/gh
```

This preparation command makes no release changes. Inspect
`publication/RELEASE_TAG/plan.json` and the `prepared/` manifests. After explicitly
authorizing replacement, repeat the command with `--publish --detach`.

Publication verifies all local checksums again before making changes. A matching
GitHub size and SHA-256 digest is reused on restart; otherwise the old asset is
deleted and its replacement uploaded under the same name. Only the exact old
asset IDs saved during preparation may be replaced. Upload responses and the
complete final inventory are checked against the prepared plan. The mosaic is
published last. API writes are sequential, with spacing and rate-limit backoff.
The uploader never changes the release object or release notes.

Progress is recorded in `publication/RELEASE_TAG/status.json` and
`publication.log`; before/after inventories are retained there. A failed upload
can temporarily leave that individual asset missing. Rerun the same publication
command to restore it and skip already verified replacements. Source files,
COGs, and previous local tile archives are retained.
