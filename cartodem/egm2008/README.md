# Reproduce the CartoDEM EGM2008 products

This workflow builds the 543 native-resolution COGs and the opaque TerrainRGB
tiles in the [CartoDEM V3R1 derived-data release](https://github.com/ramSeraph/indian_land_features/releases/tag/cartodem-30m-v3r1).
It downloads every raster input, records checksums and recipes, and can resume
completed stages. The [river method](river_flattening/README.md) describes the
experimental fitting algorithm and the measured limitations.

## Code layout

| Folder | Contents |
| --- | --- |
| This directory | Download, conversion, COG and tiling commands; dependency pins. |
| `install/` | Environment and WhiteboxTools installation scripts. |
| `river_flattening/` | River geometry, routing and fitting implementation. |
| `tests/` | Automated tests, with river tests in `tests/river_flattening/`. |
| `validation/` | Environment checks, river audits, PMTiles verification and mosaic compatibility checks. |
| `provenance/` | Recorded input/output manifests and build results. |

## Processing contract

| Stage | Result |
| --- | --- |
| Source validity | Original `-1024`, `-255`, and invalid input samples are nodata before any height calculation. |
| Datum correction | Interpret source heights as WGS84 ellipsoidal heights and compute `H = h - N` using the pinned EGM2008 geoid grid. This source-datum assumption is recorded in raster tags. |
| Ocean mask | GLO-30 WBM class 1 becomes nodata. |
| River treatment | GLO-30 WBM class 3 defines the footprint; corrected CartoDEM bank elevations supply the heights. Connected river profiles share junction elevations across sheets. |
| Other waterbodies | No new lake flattening. Class 2 is excluded from river fitting and bank samples. |
| COGs | Original grid, EPSG:4326, Float32 metres, `NoData=-255`, lossless compression and numeric overviews. |
| Display tiles | Mapbox TerrainRGB, 512 × 512, zooms 0–12, lossless WebP, quantized heights, exactly 0 m for nodata. All pixels are opaque; there is no alpha-254 marker. |
| Packaging | Per-zoom SQLite → standard MBTiles → the `pmtiles-mosaic` partition CLI. Images are copied without re-encoding. |

GLO-30 elevations are never used. Source holes on land are preserved. Stripe
correction, horizontal realignment and replacement with other DEMs are outside
this workflow.

## Software and workspace

The recorded build used Linux ARM64, Python 3.14.3, Rasterio 1.5.2 with GDAL
3.12.2/PROJ 9.8.1, WhiteboxTools 2.4.0, and the packages pinned in
[requirements.txt](requirements.txt). The display environment is pinned in
[requirements-tiles.txt](requirements-tiles.txt); the separate partitioning
environment uses [requirements-mosaic.txt](requirements-mosaic.txt).

Run commands from this directory. Choose ordinary directories on storage with
enough capacity; no machine address or mount location is required by the code.
Allow roughly 100 GiB for a fresh complete workspace, including inputs,
intermediates, COGs, encoded SQLite, MBTiles, PMTiles and temporary files. Four
CPU cores and 16 GiB RAM provide room for the four-worker tiler. River geometry
and routing use bounded, mostly sequential workers; a GPU is not used.

```bash
sudo apt-get update
sudo apt-get install -y ca-certificates curl tar gzip coreutils

dem_work=/path/to/dem-build
display_work=/path/to/display-build
tool_dir=/path/to/tools

bash install/install_pipeline_arm64.sh "$tool_dir" "$dem_work"
```

The installer creates `processing-venv`, `venv` for tiling, and `mosaic-venv`
under the selected tools directory. It downloads checksum-pinned uv and GitHub
CLI binaries, uses uv to install Python 3.14.3, installs pinned Python wheels,
and extracts checksum-pinned conda-forge WhiteboxTools/libgcc packages into the
DEM workspace. It does not install Conda or require GRASS. Package license
files and the Whitebox binary hash are retained under
`tools/whitebox-prebuilt/` in that workspace.

## 1. Download inputs and build the final COGs

```bash
"$tool_dir/processing-venv/bin/python" run_build.py \
  --workdir "$dem_work" --river-source glo30 --through cogs
```

This runs `download`, `prepare`, `geometry`, `routing`, `fit`, and `cogs` in
order. Each heavy stage runs in its own process. To inspect the stages
individually, the equivalent commands are:

```bash
"$tool_dir/processing-venv/bin/python" build.py download --workdir "$dem_work"
"$tool_dir/processing-venv/bin/python" build.py prepare --workdir "$dem_work"
"$tool_dir/processing-venv/bin/python" -m river_flattening.production geometry \
  --workdir "$dem_work" --label final --river-source glo30
"$tool_dir/processing-venv/bin/python" -m river_flattening.national routing \
  --workdir "$dem_work" --label final
"$tool_dir/processing-venv/bin/python" -m river_flattening.national fit \
  --workdir "$dem_work" --label final
"$tool_dir/processing-venv/bin/python" -m river_flattening.national apply \
  --workdir "$dem_work" --label final
```

### Download sources

[inputs.lock.json](inputs.lock.json) records URLs, byte sizes and SHA-256 hashes
for the actual inputs: 543 original CartoDEM TIFFs, 543 GLO-30 WBM TIFFs and one
geoid grid. Keep this file beside `build.py` when reproducing the build.

- CartoDEM: the [original-data release](https://github.com/ramSeraph/indian_land_features/releases/tag/cartodem-30m-v3r1-orig).
  Its checksum-pinned `listing_files.csv` supplies the original filenames and
  sizes. The 543 rasters total 14,466,784,647 bytes.
- EGM2008: [PROJ's NGA 2.5-minute geoid grid](https://cdn.proj.org/us_nga_egm08_25.tif),
  SHA-256 `4191d471eefebf24091b56dbc604353cb3b8cf8cc70e448bb9ae56a272bef17a`.
- GLO-30: water classification files from the [Copernicus DEM public S3 dataset](https://registry.opendata.aws/copernicus-dem/).
  For a sheet containing 12° N, 75° E, the URL is
  `https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N12_00_E075_00_DEM/AUXFILES/Copernicus_DSM_COG_10_N12_00_E075_00_WBM.tif`.
  `mask_name()` derives the filename from the sheet centre. The same WBM file
  supplies both ocean and river classes; no GLO-30 elevation TIFF is downloaded.

Files are saved under `inputs/cartodem/`, `inputs/geoid/`, and
`inputs/glo30-wbm/`. Existing files are checksum-verified before reuse; an
unexpected source inventory or checksum fails. An unavailable WBM is permitted
only if its entire CartoDEM sheet already has no valid elevation. All 543 masks
were available for this build.

### Datum, nodata and ocean handling

`cartodem_to_egm2008.py` masks both original fill codes first. It samples geoid
undulation bilinearly at native CartoDEM pixel centres and subtracts it from
valid elevations. The raster grid is unchanged. It samples WBM categorically
with nearest neighbour, treating class 0 as valid land even when the source
TIFF labels zero as nodata. Class 1 is removed; unknown classes and missing
mask/geoid coverage fail validation. Every removed or originally missing pixel
is stored as `-255`, with that value declared as nodata. A valid corrected
height colliding with the reserved sentinel fails explicitly.

These corrected, ocean-clipped sheets are `sheets/coastal/*.tif`. They supply
the unchanged baseline and all bank-height observations for the river stage.
The [river method](river_flattening/README.md) describes geometry, Whitebox
direction evidence, the shared-junction fit and its application.

### COG export

The `apply` stage writes the fitted values only where valid baseline pixels
fall inside supported river sections. Unresolved river sections, other
waterbodies, land and nodata retain their baseline values.
`cog.py` creates Float32 COGs with DEFLATE compression, floating-point predictor,
512-pixel blocks and averaged numeric overviews. It reopens every output and
checks every native value and validity sample against the applied raster.

Final files are `output/final/cogs/*.tif`, retaining the original sheet names.
`rivers/final/cogs.json` records file hashes, sizes, nodata and treatment counts;
`rivers/final/audit.json` records native gradients and cross-sheet checks.

## 2. Convert final COGs to TerrainRGB

```bash
"$tool_dir/venv/bin/python" rebuild_tiles.py \
  --workdir "$display_work" \
  --cog-dir "$dem_work/output/final/cogs" \
  --cog-manifest "$dem_work/rivers/final/cogs.json" \
  --workers 4 --tiles-only
```

Use `--tiles-only`: packaging is performed by the standard tool in the next
step. The tiler verifies all COG hashes, builds a shared VRT and a separate
native 0/255 validity mosaic, and resamples numeric heights onto the Web Mercator
tile grid. It averages valid samples at z0–11 and uses bilinear sampling at z12.
It never resizes packed TerrainRGB colours to construct lower zooms.

Coverage is sampled separately with nearest neighbour from the native mask,
which has no nodata metadata and no overviews. Thus averaged COG overviews
cannot silently make a sampled hole valid. A hole smaller than a coarse output
pixel can disappear when no output pixel centre samples it.

TerrainRGB decoding is `height = -10000 + 0.1 * (R*65536 + G*256 + B)`.
Valid heights use `round_bits = 16 - max(zoom, 5)`:

| Zoom | Height increment | Maximum rounding error |
| --- | --- | --- |
| 0–5 | 204.8 m | 102.4 m |
| 6 | 102.4 m | 51.2 m |
| 7 | 51.2 m | 25.6 m |
| 8 | 25.6 m | 12.8 m |
| 9 | 12.8 m | 6.4 m |
| 10 | 6.4 m | 3.2 m |
| 11 | 3.2 m | 1.6 m |
| 12 | 1.6 m | 0.8 m |

Nodata bypasses quantization and is set to RGB `(1,134,160)`, exactly 0 m.
The internal RGBA array is fully opaque; the inspected lossless WebP outputs
decode as RGB with no used alpha channel. No 254 or 0.1 m sentinel is present.
Consumers use ordinary Mapbox TerrainRGB decoding with `tileSize: 512`.
Use the COGs when actual nodata distinction or unquantized elevations matter.

The encoded images are stored in resumable `tiles/display/z0.sqlite` through
`z12.sqlite`. Every image is decoded and compared with its expected bytes;
valid heights are checked against the quantization error bound. The four-worker
display build produced 227,426 tiles in about 2 h 15 min, including its original
packaging stage. See [TILING_SETUP.md](TILING_SETUP.md) for tiling details.

## 3. Export MBTiles and partition with pmtiles-mosaic

```bash
"$tool_dir/venv/bin/python" export_mbtiles.py --workdir "$display_work" \
  --outdir "$display_work/publication/mbtiles"
bash partition_mbtiles.sh "$display_work/publication/mbtiles" \
  "$display_work/publication/pmtiles" "$tool_dir/mosaic-venv/bin/partition"
"$tool_dir/mosaic-venv/bin/python" -m validation.verify_partitioned \
  --workdir "$display_work" --pmtiles-dir "$display_work/publication/pmtiles"
```

Export copies the existing WebP blobs into the standard MBTiles schema and
converts XYZ `y` to TMS `tile_row = 2**z - 1 - y`. All image hashes are checked.
`partition_mbtiles.sh` calls **pmtiles-mosaic 0.2.0**, supplying all thirteen
MBTiles files with repeated `--from-source`, a 2,000,000,000-byte limit, and
`--no-cache`. The tool creates geographic partitions and its own mosaic JSON.
The final verifier checks every packaged tile against the original databases,
including missing, duplicate and changed tiles.

This build produced four archives named
`Bhuvan_CartoDEM_v3r1_TerrainRGB-part0000.pmtiles` through `part0003.pmtiles`,
totaling 7,273,487,442 bytes. The IndianOpenMaps production partition-aware
resolver was also checked against every tile. To repeat that application check:

```bash
node validation/verify_mosaic.mjs /path/to/pmtiles-mosaic/dist \
  "$display_work/publication/pmtiles/Bhuvan_CartoDEM_v3r1_TerrainRGB.mosaic.json" \
  "$display_work/publication/pmtiles/.validation/partition-info.json"
```

The private index used here is generated by `validation/verify_partitioned.py` for this
validation and is not uploaded. See the publication commands in `TILING_SETUP.md`.
The PMTiles writer formerly embedded in `production_tiles.py` is a historical
code path and is not used by these commands.

## Resume, validation and publication

`run_build.py --resume-from STAGE --through cogs` resumes a DEM stage. Adding
`--detach` writes a durable log and status record; wait for completion before
starting the display build. A bounded river pilot can use
`--label pilot-final --sheets cdnd43v.tif cdnd43w.tif --through cogs` in a fresh
workspace. The download stage still verifies the complete original inventory.

Geometry, routing evidence and fitted components are checkpointed. Input and
recipe hashes reject incompatible caches. Temporary routing rasters are removed
after success; their hashes, timing and aggregated direction evidence remain.
Rebuilding them requires rerunning routing, not filling the published COGs.
The display tiler resumes completed x-column batches. MBTiles export verifies
completed files before reuse; an incomplete `.mbtiles.part` must be inspected
and removed before retrying that zoom.

Run the synthetic and small-raster tests in the processing environment:

```bash
"$tool_dir/processing-venv/bin/python" -m unittest discover -s tests -t . -v
```

Tests cover both nodata codes, geoid correction, ocean masking, native validity,
COG preservation, direction conventions, shared junctions, connected fitting,
cross-sheet application, tile resampling, quantization, MBTiles row conversion,
and byte preservation. They do not establish
that every experimental river surface is correct; the national audit's measured
limitations are documented in the river method.

The rebuild commands do not publish. [TILING_SETUP.md](TILING_SETUP.md) uses
`gh-release-tools` for data uploads and `listing_files.csv`. Publish only the
COGs, PMTiles, mosaic JSON and listing; build records stay in the repository or
workspace. The recorded COG manifest and compact results are in
[provenance/](provenance/README.md). Runtime provenance metadata can change a
rerun's file hashes; input checksums, native pixels and the processing settings
are the basis for reproducing the data.
The [release notes](../RELEASE_NOTES.md) describe the published products.
Downloaded rasters, working databases, credentials and machine connection
settings are excluded from the commit.
