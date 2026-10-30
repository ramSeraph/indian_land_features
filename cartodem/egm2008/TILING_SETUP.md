# Rebuild display tiles and publish the data

This is the shorter path when the final EGM2008 COGs already exist. For downloads,
ocean masking and experimental river fitting, use [README.md](README.md).
The recorded COG manifest is [provenance/cogs.json](provenance/cogs.json).

## Install and rebuild

Run from this directory on Linux ARM64:

```bash
tool_dir=/path/to/tools
display_work=/path/to/display-build
cog_dir=/path/to/final/cogs
cog_manifest=/path/to/cogs.json

bash install_tiles_arm64.sh "$tool_dir"

"$tool_dir/venv/bin/python" rebuild_tiles.py \
  --workdir "$display_work" --cog-dir "$cog_dir" \
  --cog-manifest "$cog_manifest" --workers 4 --tiles-only
```

The installer uses checksum-pinned uv and GitHub CLI binaries, Python 3.14.3,
[requirements-tiles.txt](requirements-tiles.txt), and a separate
[pmtiles-mosaic environment](requirements-mosaic.txt). Rasterio supplies GDAL
and PROJ; Pillow supplies lossless WebP. No rio-rgbify or GDAL command-line
installation is needed.

Use a fresh display workspace when changing the input COGs or encoding recipe.
The runner checks all 543 input COG hashes and does not modify them. Four CPU
cores and 16 GiB RAM provide room for four workers; one worker uses one CPU.
Reserve at least 30 GiB beyond the COGs for encoded SQLite, standard MBTiles,
PMTiles, validity masks and temporary files. The recorded four-worker build
completed in about 2 h 15 min, including its original packaging stage.

Images are 512 × 512 lossless WebP at z0–12. Valid numeric heights are averaged
at z0–11 and sampled bilinearly at z12, then quantized with
`round_bits = 16 - max(z, 5)`. At z12 that is a 1.6 m increment and at most
0.8 m rounding error. Missing samples bypass quantization and become
RGB `(1,134,160)`, exactly 0 m. All pixels are opaque. The WebP encoder omits
unused transparency; no alpha-254 sentinel or special browser decode protocol
is required. True zero and nodata are indistinguishable in the display tiles.

A separate native validity raster controls coverage by nearest pixel centre,
independently of averaged COG overviews. This does not guarantee that every
subpixel hole remains visible at coarse zooms. COGs retain normal `NoData=-255`.

The tiler resumes completed x-column batches. Zooms 0–4 run serially; z5–12 use
isolated workers with one GDAL thread each and bounded row batches. Changing
worker count does not change the encoding recipe.

## Standard MBTiles and PMTiles

```bash
mbtiles_dir="$display_work/publication/mbtiles"
pmtiles_dir="$display_work/publication/pmtiles"

"$tool_dir/venv/bin/python" export_mbtiles.py \
  --workdir "$display_work" --outdir "$mbtiles_dir"

bash partition_mbtiles.sh "$mbtiles_dir" "$pmtiles_dir" \
  "$tool_dir/mosaic-venv/bin/partition"

"$tool_dir/mosaic-venv/bin/python" -m validation.verify_partitioned \
  --workdir "$display_work" --pmtiles-dir "$pmtiles_dir"
```

Export changes the database schema and converts XYZ rows to MBTiles TMS rows.
It copies image bytes and verifies their SHA-256 values. The standard
`pmtiles-mosaic==0.2.0` CLI reads all thirteen MBTiles sources through repeated
`--from-source` arguments and chooses geographic partitions with a
2,000,000,000-byte limit. It uses `--no-cache`; temporary files stay in the
selected output directory. No image encoding or raster resampling is repeated.

The verifier checks every tile against the original encoded databases. It writes
internal build manifests and a private `.validation/partition-info.json`.
Check application compatibility with the same JavaScript mosaic resolver used
by the serving application:

```bash
node validation/verify_mosaic.mjs /path/to/pmtiles-mosaic/dist \
  "$pmtiles_dir/Bhuvan_CartoDEM_v3r1_TerrainRGB.mosaic.json" \
  "$pmtiles_dir/.validation/partition-info.json"
```

For this build, all 227,426 tiles were checked with IndianOpenMaps' production
partition-aware resolver. Four PMTiles archives total 7,273,487,442 bytes and
retain the `Bhuvan_CartoDEM_v3r1_TerrainRGB-part0000.pmtiles` through
`part0003.pmtiles` filenames.

## Publish with gh-release-tools

Publish the COGs, PMTiles, the tool-generated mosaic JSON and the TIFF listing.
Build manifests, checksum files, software records, source archives and the private
partition index are not release assets. Keep them in the workspace or repository.

Authenticate GitHub CLI, choose the intended existing release, then use
[gh-release-tools](https://github.com/ramSeraph/release_tools). These commands
replace matching filenames:

```bash
export PATH="$tool_dir/bin:$PATH"
repo=OWNER/REPOSITORY
release_tag=RELEASE_TAG

gh auth status

uv tool run --from gh_release_tools==0.1.14 upload-to-release \
  --repo "$repo" --release "$release_tag" --folder "$cog_dir" \
  --extension .tif --batch-size 1 --overwrite

uv tool run --from gh_release_tools==0.1.14 upload-to-release \
  --repo "$repo" --release "$release_tag" --folder "$pmtiles_dir" \
  --extension .pmtiles --batch-size 1 --overwrite

gh release upload "$release_tag" \
  "$pmtiles_dir/Bhuvan_CartoDEM_v3r1_TerrainRGB.mosaic.json" \
  --repo "$repo" --clobber

uv tool run --from gh_release_tools==0.1.14 generate-lists \
  --repo "$repo" --release "$release_tag" --extension .tif
```

`generate-lists` reads the actual release inventory and replaces
`listing_files.csv` with its name, size and download-URL columns. The finished
release contains 543 COGs, four PMTiles files, one mosaic JSON and one listing:
549 assets. There is no public partition-info file.

Release notes should identify the CartoDEM source, EGM2008 correction,
GLO-30 ocean masking and **experimental river flattening using the GLO-30 river
mask**, and link to the reproduction commit. The release text is in
[../RELEASE_NOTES.md](../RELEASE_NOTES.md). Data uploads and listing generation
do not edit the release notes.
