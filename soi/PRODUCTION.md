# SOI contour-derived elevation sheets — WIP

This workflow produces intermediate 30 m elevation sheets from the released
Survey of India contour vectors. The sheets will be replaced by final data.
**Waterbody and river flattening are deferred and have not been applied.**

## Processing

1. Download the exact assets in `manifests/soi-contours.json`, resuming incomplete
   transfers and checking byte counts and SHA256 digests.
2. Use the WGS84 polygons in `index_50k.geojson` and the selected sheets in
   `manifests/sheets_interior.txt`. Extract complete intersecting features from
   all matching local Parquet partitions with iomaps, retaining their `VALUE`
   heights and uncut geometry for classification.
3. Identify circular and triangular height symbols, infer the contour interval
   separately for each sheet, exclude supported form-line candidates and apply
   bounded height recovery. Supported suspect heights become unknown constraints
   for interpolation; unresolved contradictory constraints stop that sheet.
4. Interpolate in the sheet's WGS84 UTM zone using the Python Imagico adaptation.
   **All accepted height points are included.** Production skips the holdout
   experiment; calibration and comparison scripts are separate.
5. Build neighbouring first passes on demand. Solve shared edge patches using a
   2 km replacement band on each side and 3.6 km of additional context. Solve
   intersecting edge strips at corners after applying their shared edge patches.
6. Validate the sheet mask, grid, all fixed constraint pixels and all height
   points. Export a losslessly compressed Cloud Optimized GeoTIFF and verify
   that its full-resolution pixels exactly equal the checked raster.
7. Retain the classified circle/triangle points in `data/height_points/` as one
   GeoJSONL per sheet. Retain compact local reports and logs, release temporary
   job files, and upload only the elevation COG when a release is specified.

No scanned toposheet, GLO-30 elevation or manual review file is a production
elevation input. WGS84 is the agreed horizontal grid. Source heights are retained
without an assumed vertical-datum conversion; these outputs are not labelled as
a verified EGM2008 conversion.

## Reproduce the inputs

Run commands from `soi/`. Use Python 3.12.9, `uv`, `curl`, and GDAL suitable for
the locked Rasterio build. The tested GDAL version is 3.11.4. Python dependencies
are pinned in `produce_sheets.py.lock`.

```sh
python3 download_inputs.py manifests/soi-contours.json data/contours/soi-raw
```

The snapshot contains four Parquets, 4.518 GiB and 7,362,891 rows. Downloaded files
are renamed from `.part` only after verification. The metadata asset is pinned
as well. Source geometries and heights are not edited during download.

The production tool needs the local-Parquet extraction option added to iomaps.
The base commit and repository are recorded in `vendor/iomaps-source.json`, with
the exact changes in `vendor/iomaps-local-input.patch`:

```sh
git clone https://github.com/ramSeraph/indianopenmaps.git /path/to/indianopenmaps
git -C /path/to/indianopenmaps checkout 20b8154a0bdeec49fd1db5a5986c436faccd3c93
git -C /path/to/indianopenmaps apply /path/to/soi/vendor/iomaps-local-input.patch
uv sync --frozen --project /path/to/indianopenmaps/python
```

Install the DuckDB spatial extension in that environment. Local extraction does
not need the remote catalogue or HTTP extension. The run records installed
versions and hashes of the iomaps extraction and filtering modules.

The inventory can be regenerated from the local Parquets:

```sh
uv run inventory_contour_sheets.py \
  data/contours/soi-raw/SOI_Contours_raw.parquet.meta.json \
  --output-dir data/contours/soi-raw/sheet-inventory --plot
```

Use a new output directory. The inventory intersects actual geometries with the
5,113 WGS84 sheet polygons: 2,321 sheets have linework beyond a 2 m inward edge
band, 83 have edge-only linework, 21 have point contacts only and 2,688 have none.
Presence does not establish complete coverage or consistent elevations.

## Run and resume

For a limited local run using the tested three-sheet selection:

```sh
uv run --locked produce_sheets.py \
  --sheets-file calibration/production_overlap_sheets.txt \
  --iomaps ../../indianopenmaps/python/.venv/bin/iomaps \
  --output-dir data/contour_dem/local-check --limit 3
```

For the full selection and an existing release:

```sh
uv run --locked produce_sheets.py \
  --sheets-file manifests/sheets_interior.txt \
  --iomaps ../../indianopenmaps/python/.venv/bin/iomaps \
  --output-dir data/contour_dem/production \
  --release soi-dem-30m --prune-uploaded
```

Omit `--release` and `--prune-uploaded` to keep outputs locally. `--limit N`
processes at most N additional sheets; rerun the same command without the limit
to continue. `--priority 56M_9 83E_11` changes queue order within the selected
sheets. Create `OUTPUT/STOP` to stop after the current sheet; remove it to resume.
`--no-overlap` requests independent first-pass outputs explicitly.

Processing is sequential, with a 3 GiB free-space reserve, 30-minute sheet
timeout, and bounded caches of 16 first passes and 64 overlap patches. An
exclusive run lock prevents concurrent writers. Changes to inputs, code,
software versions or numerical options require a new output directory.
Failed sheets are recorded locally and produce no release asset.

Uploads use gh-release-tools 0.1.14. Completed outputs are removed locally only
after remote size and SHA256 verification. Existing remote names are not
overwritten. `run_production_batch.py` wraps the producer and refreshes the
release listing after the selected queue finishes.

## Outputs and remaining limitations

Release assets are `SOI_DEM_30m_SHEET.tif`. Code, dependency locks, documentation,
source manifests and calibration scripts are committed in this repository.
There are no workflow, validation or per-sheet diagnostic ZIPs in the release.

Each COG uses the sheet's WGS84 UTM zone, Float32 elevations in metres, DEFLATE
compression, 512-pixel tiles and averaged overviews. Local `run.json` records the
input and implementation identities; `state.json`, `reports/` and `logs/` support
resuming and diagnosing the queue. The separate point GeoJSONL retains source
attributes and circle/triangle classifications. Triangles are triangulated
height stations, not a general peak classification.

Shared patches preserve measured heights and leave pixels outside their
replacement region unchanged. A two-pixel outer frame uses the first-pass
surface as a model boundary condition. Both sheet owners apply the same saved
patch. Corner patches cover the actual intersection of incident edge strips,
including the rotation between geographic boundaries and UTM grids.

The COG tag `SOI_SEAM_STATUS` distinguishes complete, partial and isolated
correction. Missing neighbours, source conflicts, cross-UTM joins and unresolved
unknown-contour identities remain explicit limitations. Unknown groups crossing
a patch boundary are rejected. The accepted local correction is not claimed to
equal a single country-wide interpolation or to bound all neighbouring influence.

## Local validation

The five successful regional pilots included all 589 classified height points
and preserved 894,953 fixed constraint pixels. Comparison against Standard
Copernicus GLO-30, before overlap correction, gave:

| Sheet | Height points | Mean SOI minus GLO-30 | RMSE |
| --- | ---: | ---: | ---: |
| 56M_9 | 104 | -3.08 m | 6.44 m |
| 56M_10 | 115 | -1.73 m | 5.52 m |
| 45F_4 | 150 | -2.97 m | 4.92 m |
| 73D_15 | 122 | -1.18 m | 7.02 m |
| 83E_11 | 98 | -13.72 m | 23.65 m |

These are differences from another DEM, not independent survey accuracy.
53A_16 remained a failed case with 14 conflicting sample locations; its source
checks were not bypassed. Other incomplete or inconsistent sheets can fail too.

The integrated 56M_9/56M_10/56M_14 overlap rehearsal preserved all 315 height
points and 143,441 fixed pixels. Within 300 m of the joins, GLO-30 RMSE improved
from 5.48 m to 4.64 m. Saved patch inputs replayed to identical corrections.
Output COGs were 2.02–2.13 MiB. With the full selection, a real four-sheet corner
between 55P_8, 55P_12, 56M_5 and 56M_9 also passed, preserving its 503 fixed
pixels. Sheet 56M_9 remained marked partial because other neighbours failed.

Calibration commands, historical decisions and detailed observations remain
separate in [calibration/NOTES.md](calibration/NOTES.md) and
[calibration/PRODUCTION_NOTES.md](calibration/PRODUCTION_NOTES.md). The calibration
scripts have their own dependency locks. GLO-30 is used only for comparison.

The earlier checks passed 124 DEM/comparison tests, three inventory tests and
14 iomaps extraction tests, followed by four numerical patch tests and five
overlap integration tests. The inventory tests run separately from Rasterio
tests because loading both native library bundles together caused duplicate
Arrow filesystem registration in the tested environment.

## Attribution

The Imagico numerical port and raster adapter are GPL-3.0-or-later, matching
Christoph Hormann's original code. See [LICENSE.imagico.txt](LICENSE.imagico.txt),
[imagico_provenance.json](imagico_provenance.json) and `contour_interpolate-0.3/`.
NumPy/SciPy filtering, missing-data masks and constraint protection mean this is
not a bit-for-bit CImg reproduction. The source data's licensing statement is
retained in [DATA_LICENSE.source.md](DATA_LICENSE.source.md); credit Survey of
India and the Datameet community for the source data.
