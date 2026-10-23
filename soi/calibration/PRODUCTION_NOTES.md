# Historical production experiments

These notes record the earlier packaging and validation experiments. For the current workflow, see [PRODUCTION.md](../PRODUCTION.md). Commands assume `soi/` as the working directory.

# SOI contour-derived elevation sheets

This workflow produces one 30 m elevation sheet from each selected SOI 1:50,000
sheet with interior contour data. The source snapshot has 2,321 candidate sheets.
Presence is not a guarantee of complete coverage or consistent elevations. A
failed sheet keeps its source extract, classified points and failure diagnostics;
it does not produce an elevation GeoTIFF.

## Production stages

1. Download the exact assets in `manifests/soi-contours.json`, resuming incomplete
   transfers and verifying each byte count and SHA256.
2. Use the WGS84 polygons in `index_50k.geojson` and the pinned selection in
   `manifests/sheets_interior.txt`. `inventory_contour_sheets.py` can regenerate
   the inventory from the local Parquets, using actual line intersections and a
   2 m inward edge band.
3. Extract complete intersecting features from all matching local partitions
   with `iomaps cli extract --parquet`. Preserve original `VALUE` attributes and
   the uncut geometries for classification and later boundary processing.
4. Apply the existing circle/triangle classifiers, per-sheet interval heuristic,
   form-line exclusion, bounded height-recovery rules and contact handling.
   Supported suspect heights become unknown constraints for Imagico; unresolved
   contradictory heights stop that sheet. No universal contour interval is used.
5. Interpolate in the sheet's WGS84 UTM zone with the Python Imagico adaptation.
   **All available height points are included in this one production pass.**
   Production uses `--skip-holdout`; held-out-point calibration is separate.
6. Build neighbouring first passes on demand and reconcile shared edges with
   one saved patch per edge. Reconcile intersecting edge strips at corners after
   applying their shared edge patches. Each owner uses the same saved pixels;
   processing order does not create independently estimated copies of a seam.
7. Check the exported mask, CRS, resolution, valid-cell count, every fixed
   constraint pixel and every height point. Export a losslessly compressed COG
   and check that its full-resolution pixels exactly equal the original raster.
8. Package the source extract, classified points, original/inferred contour
   records, interpolation grids, diagnostic rasters, logs and provenance into a
   separate per-sheet audit ZIP. Include the original DEM and the complete
   numerical inputs, outputs and provenance of applied overlap jobs. Check its
   CRC before removing staging files.
9. Upload the GeoTIFF and audit ZIP through **gh-release-tools 0.1.14**, with
   `--create-extra-releases`. Verify GitHub's asset size and SHA256 before pruning
   this run's local packages. Existing remote names are never overwritten.

No scanned toposheet, GLO-30 elevation or manual review file is a production
elevation input. WGS84 is the agreed horizontal grid; the original source heights
are preserved without an assumed vertical-datum correction. The result is not
labelled as a verified EGM2008 conversion.

## Reproduce the downloads

Use Python 3.12 or later, `uv`, `curl`, and an installed GDAL suitable for the
locked Rasterio build. The tested GDAL is 3.11.4. Exact Python packages are locked
in `produce_sheets.py.lock`; the run manifest also records their installed
versions, Python build, GDAL version and implementation hashes.

```sh
python3 download_inputs.py manifests/soi-contours.json data/contours/soi-raw
```

`inputs_verified.json` records the verified inputs, URLs, checksums and local
file timestamps. A cached file with the wrong size or checksum is rejected.
Downloaded files are atomically renamed from `.part` only after verification.

The source snapshot occupies 4.518 GiB and contains 7,362,891 rows in four
Parquets. The metadata file is also pinned. No reprojection or editing is applied
to the downloaded sources.

## Reproduce the iomaps extraction tool

The production command needs the local-Parquet option added to iomaps. Its base
Git commit and repository are recorded in `vendor/iomaps-source.json`; the exact
change is `vendor/iomaps-local-input.patch`. In a fresh checkout:

```sh
git clone https://github.com/ramSeraph/indianopenmaps.git /path/to/indianopenmaps
git -C /path/to/indianopenmaps checkout 20b8154a0bdeec49fd1db5a5986c436faccd3c93
git -C /path/to/indianopenmaps apply /path/to/soi/vendor/iomaps-local-input.patch
uv sync --frozen --project /path/to/indianopenmaps/python
```

The DuckDB spatial extension must be installed for that environment. Local
extraction skips the remote catalogue and HTTP extension. The producer records
the installed iomaps version and hashes of its extraction/filter modules.

## Run and resume

The main release must exist before using `--release`. Release creation is a
separate action performed only after the limited local checks are reviewed.

```sh
uv run --locked produce_sheets.py \
  --output-dir data/contour_dem/production-v1 \
  --sheets-file manifests/sheets_interior.txt \
  --iomaps ../../indianopenmaps/python/.venv/bin/iomaps \
  --release soi-dem-30m --prune-uploaded
```

Omit `--release` and `--prune-uploaded` for local output. `--limit 2` stops cleanly
after two additional sheets; rerun without the limit to continue the same run.
`--priority 56M_9 83E_11` changes queue order without changing the selection.
`--no-overlap` explicitly requests independent first-pass output. The default
includes the accepted 2 km overlap correction, with 3.6 km of additional context.

The default is **one sheet at a time** with a **3 GiB free-space reserve** and a
30-minute timeout per sheet. The runner holds an exclusive lock, checks for
changed source files and code between sheets, and checkpoints after packaging
and verified upload. Reusing a run with different source identities, code,
software versions, index or options is refused. A changed build needs a new run
directory. Interrupted work is retained under `interrupted/` before retrying;
the runner refuses to race a still-running previous sheet process.

Create `OUTPUT/STOP` to stop after the current sheet. Remove that marker before
resuming. The runner also stops cleanly when the disk reserve is reached.
`run_production_batch.py` runs this same producer with verified uploads and
pruning, then calls gh-release-tools `generate-lists` only after the queue has
finished. It accepts the same output folder, selection, repository, release and
priority sheets. It does not run water conditioning.
It processes nearby sheets in 0.5-degree blocks and retains at most 16 completed
first passes and 64 completed overlap jobs in reusable caches. The current
neighbourhood can temporarily exceed those limits while a sheet is being built.
Unused cache entries can be discarded and deterministically regenerated; source
Parquets and verified release packages are not treated as caches.

Files in the output directory:

| Path | Contents |
| --- | --- |
| `run.json` | Pinned configuration, code/environment hashes and source verification |
| `state.json` | Queue progress and per-sheet outcome |
| `reports/SHEET.json` | Local checks or failure details |
| `logs/SHEET.log` | Extraction and interpolation log |
| `initial-cache/`, `patch-cache/` | Bounded reusable first passes and canonical edge/corner patches |
| `initial-results/` | First-pass outcomes, including failed neighbours |
| `receipts/SHEET.json` | Verified remote URLs, sizes, SHA256 and asset IDs |
| `packages/SHEET/SOI_DEM_30m_SHEET.tif` | Elevation COG pending upload/pruning |
| `packages/SHEET/SOI_DEM_30m_SHEET.audit.zip` | Source constraints and processing evidence |
| `packages/SHEET/SOI_DEM_30m_SHEET.failed.zip` | Evidence for a sheet which could not be built |

To regenerate a sheet independently, extract its audit ZIP and run
`contours_to_dem.py` against `source/SHEET/contours.geojsonl`, using the pinned
index and `--skip-holdout`. Paths and elapsed timings in audit JSON reflect the
original run and need not be byte-identical on another machine. Numerical
reproduction should be checked on the full-resolution raster with the recorded
software environment.

## Calibration is separate

`calibration/run_dem_pilot.py` builds a limited set with held-out-point evaluation
and then uses `compare_dem.py` to compare each final all-point DEM with standard
Copernicus GLO-30. This is **not part of the production queue**. Existing
`audit_*.py` scripts and `test_*.py` files are retained as diagnostics/tests.

```sh
uv run calibration/run_dem_pilot.py \
  --sheets 56M_9 45F_4 73D_15 83E_11 56M_10
```

GLO-30 windows preserve their native pixel grid and are then resampled bilinearly
onto the exact SOI output grid. Reports record URLs, native grids, cached-window
checksums, parameters and the comparison-script hash. No fitted horizontal or
vertical correction is applied. Differences are **SOI minus GLO-30**. GLO-30 is a
surface model, so vegetation, buildings, epoch, registration and interpolation
all contribute to differences; these are not independent survey error estimates.

## Local validation on 2026-10-23

| Sheet | Supported sheet cells | Height points included | Mean difference | RMSE |
| --- | ---: | ---: | ---: | ---: |
| 56M_9 | 100% | 104 | -3.08 m | 6.44 m |
| 56M_10 | 100% | 115 | -1.73 m | 5.52 m |
| 45F_4 | 97.57% | 150 | -2.97 m | 4.92 m |
| 73D_15 | 100% | 122 | -1.18 m | 7.02 m |
| 83E_11 | 100% | 98 | -13.72 m | 23.65 m |

All **589 height points** and **894,953 fixed raster constraint pixels** survive
the final output at Float32 precision. The four sheets with earlier remote
extracts match those source features, geometries and attributes exactly; 56M_10
is a new extract. The initial boundary validator incorrectly projected only
polygon corners. Densifying the geographic sheet edges before UTM projection
fixed that check; no DEM was changed to make it pass.

53A_16 remains a deliberately included difficult case: it stops at **14
conflicting sample locations**, preserving its classified points and diagnostics.
This has not been bypassed or labelled a successful DEM. The automated batch
must be expected to report other source conflicts and incomplete coverage.

The production-only rehearsal on 56M_9 and 83E_11 produced pixel-identical rasters
to their all-point pilot outputs, without rerunning the holdout stage. Build/check
times were 34.5 and 42.9 seconds; combined COG/audit sizes were approximately
5.55 and 19.45 MiB. These are individual observations, not a country-wide runtime
or storage guarantee.

Machine-readable pilot results are in
`data/contour_dem/production-pilot/pilot_status.json`, with per-sheet comparisons
under `data/comparisons/production-pilot/`. The 124 DEM/comparison regression
tests, three Parquet inventory tests and 14 iomaps extraction tests passed.
Arrow inventory tests run separately from the GDAL/Rasterio suite because loading
both packaged native libraries in one process caused duplicate Arrow filesystem
registration in this environment.

## Boundary and water stages

The initial interpolation is an independent-sheet, unconditioned DEM. The
default producer now applies local overlap corrections before packaging the
final COG. Audit packages retain the first pass as `dem/dem_initial.tif`, its
report, source features and constraints, plus `overlaps/JOB/inputs.npz`,
`patch.npz` and `report.json` for every applied numerical job.

The output TIFF tag `SOI_SEAM_STATUS` and `dem/overlap.json` distinguish complete
correction for available selected neighbours, partial correction and an isolated
sheet. `seam_reconciled` stays false when a required job is blocked or a selected
neighbour failed. The report names the unresolved edges/corners. A cross-UTM
join or unknown contour crossing a patch boundary is flagged, not silently
resampled or relabelled corrected. Missing source coverage remains missing.

**Waterbody and river flattening are deferred at the user's request while a
separate agent researches them.** They are not applied by this workflow. The
following previously recorded options and downloaded masks are retained as
background for that later work; no flattening experiment is currently running.

Water conditioning would be a separate derived output, preserving the raw DEM:

- CARTODEM v3r1 waterbody polygons identify candidate lake footprints. Water
  heights must be derived from SOI elevations, with unsupported shorelines,
  conflicting height constraints and cross-sheet features explicitly handled.
- GLO-30 WBM **class 3** identifies rivers. River cross-sections should be level
  while longitudinal profiles preserve downstream fall and junction continuity.
  A whole connected river network must not be assigned one constant elevation.
- Categorical masks use nearest-neighbour alignment. WBM class **0 is valid
  land**, even when tagged as nodata; unavailable mask coverage is a separate
  state. GLO-30 elevation values must not be silently substituted into the SOI DEM.
- Compare raw and conditioned variants separately against GLO-30. A comparison
  using GLO-30-derived masks is not entirely independent of that reference.

The waterbody Parquet is pinned in `manifests/cartodem-waterbodies.json`. Nine
pilot WBM tiles are pinned in `manifests/glo30-wbm-pilot.json`; they are available
under each AWS tile's `AUXFILES/` path. They are downloaded separately from the
contour production pass:

```sh
python3 download_inputs.py manifests/cartodem-waterbodies.json data/hydro/cartodem
python3 download_inputs.py manifests/glo30-wbm-pilot.json data/hydro/glo30-wbm
```

River and waterbody conditioning remain separate from production. Publishing
these unconditioned sheets does not claim either stage is complete.

### Local overlap implementation and checks

`seam_interpolate.py` now implements a local shared-constraint solve. The
calibration adapter `calibration/reconcile_adjacent_pilot.py` assembles retained
constraint grids from neighbouring sheets on their common UTM pixel lattice.
It runs Imagico across the join and replaces a 2 km band on either side, with
3.6 km of additional constraint context. A two-pixel outer frame retains the
first-pass elevations as explicitly labelled model boundary conditions. Every
measured height takes priority. The corrected sheets are exported as COGs,
checked for pixel equality, and kept separately from the original pilot DEMs.

The user accepted this local approach and its tested 2 km width for production.
The replacement width remains **not a universal bound** on
the influence of neighbouring information. This is not a claim that local
correction equals a single country-wide interpolation. Unknown-height groups
inside the patch are re-estimated; a group extending outside the patch requires
an expanded region. Unknown groups meeting the sheet join are rejected until
their original vector identities can be reconciled. Different UTM zones are
also rejected by this pilot adapter rather than silently resampling constraints.

Reproduce the two-sheet test after building the original pilot sheets:

```sh
uv run --locked calibration/reconcile_adjacent_pilot.py \
  data/contour_dem/production-pilot/56M_9 \
  data/contour_dem/production-pilot/56M_10 \
  --output-dir data/contour_dem/overlap-reproduction
uv run --locked calibration/evaluate_seam_pilot.py \
  data/contour_dem/overlap-reproduction
```

The reference comparison uses the cached, aligned GLO-30 windows from the
original pilot; GLO-30 heights never enter the overlap solve. Input hashes,
implementation hashes, package versions, solver options, inferred heights and
COG/constraint checks are recorded in the patch report. Calibration code and its
dependency locks remain separate from the production queue.

On 56M_9/56M_10, the local solve plus export/checks took **12.4 seconds** and
preserved all **94,284 measured constraint pixels** across the two outputs.
Pixels outside the replacement strip are unchanged. Within 300 m of the join,
RMSE against GLO-30 improved from **6.15 m to 5.49 m**. The RMS difference between
the two DEMs' adjacent-pixel steps across the join improved from **5.06 m to
0.84 m**; subtracting GLO-30's local step avoids treating real terrain slope as a
zero-height seam. These are comparison metrics, not independent survey accuracy.
The 1,761 tested pixel steps across the outer replacement boundary were exactly
unchanged, so no new jump was introduced where the patch meets the old surface.

Results, corrected COGs and the comparison image are under
`data/contour_dem/seam-pilot-56M_9-56M_10-v2/`. Four regression tests cover
preserved heights, unchanged surroundings, unknown-group handling and a
controlled four-quadrant corner. The first two proposed real four-sheet groups
could not be completed: 56M_6 stopped at two conflicting sample locations
(510/517 m), and 56M_13 stopped at one. Their failed builds and diagnostics remain in
`data/contour_dem/production-pilot/`; neither was bypassed to complete a corner.

The three available sheets 56M_9, 56M_10 and 56M_14 were also solved together,
including horizontal and vertical joins meeting at the missing 56M_13 corner.
This took **45.3 seconds**, preserved **143,441 measured constraint pixels** and
left the unavailable sheet masked. Within 300 m of these joins, RMSE against
GLO-30 improved from **5.48 m to 4.62 m**; adjacent-pixel step residual RMS fell
from **5.29 m to 0.65 m**. All 3,484 outer replacement-boundary steps were
unchanged. Results are in `data/contour_dem/seam-pilot-three-sheets/`. Reproduce
with the same commands above, adding the 56M_14 pilot folder to the interpolation
command and choosing a fresh output directory. This tests an incomplete
neighbourhood; it does not substitute for the complete four-sheet-corner check.

`sheet_overlaps.py` and `overlap_cache.py` connect the solver to
`produce_sheets.py`. Corner patches cover the actual intersection of incident
edge strips, including the rotation of UTM relative to geographic sheet edges.
The numerical input archives allow a patch to be recomputed with
`seam_interpolate.reconcile_patch` without downloading neighbouring outputs.
Cross-UTM joins and unresolved unknown-contour identities remain explicitly
flagged limitations. Five production integration tests cover canonical corner
jobs, processing-order independence, unavailable neighbours, cache reuse and
cross-UTM rejection and low-disk cache eviction, in addition to the four numerical
patch tests.

The local integration rehearsal uses:

```sh
uv run --locked produce_sheets.py \
  --sheets-file calibration/production_overlap_sheets.txt \
  --output-dir data/contour_dem/production-overlap-rehearsal --limit 3
```

This limited selection contains 56M_9, 56M_10 and 56M_14. The final country-wide
run uses `manifests/sheets_interior.txt`, so a failed selected neighbour such as
56M_13 is recorded as unavailable rather than omitted from its expected jobs.

This integrated rehearsal passed for all three sheets, preserving **315 height
points** and **143,441 fixed constraint pixels**. All three archived patch inputs
replayed to pixel-identical corrections. Within 300 m of the joins, RMSE against
GLO-30 improved from **5.48 m to 4.64 m**. COG sizes were 2.02–2.13 MiB and audit
archives 8.45–9.27 MiB. The first output took 140 seconds including all three
first-pass builds and shared patches; the remaining sheets reused those caches.
Resuming the completed run skipped all three outputs successfully. The subsequent
disk-pressure safeguard was covered by the integration test; the numerical
solver and rasterization were unchanged.

The first run using the full 2,321-sheet selection subsequently solved a complete
real corner between **55P_8, 55P_12, 56M_5 and 56M_9**, preserving its 503 fixed
constraint pixels. Sheet 56M_9 passed all checks, including its 104 height points
and 43,826 fixed constraint pixels. It is correctly marked `partial` because
other expected neighbours (55P_16, 56M_6 and 56M_13) failed their source checks.
The complete corner and the available three-sheet corners were applied; the
unavailable joins were retained as explicit gaps in the correction report.
All nine numerical patch archives for this output replayed exactly. Across its
four expected edge bands (within 300 m, including the unavailable eastern
neighbour), RMSE against GLO-30 changed from **6.39 m to 6.11 m**.

### Release archives

`package_workflow.py` creates `SOI_DEM_30m_workflow.zip` with the production code,
dependency lock, pinned source manifests, sheet grid/selection, iomaps patch,
licenses and this document. `SOI_DEM_30m_validation.zip` keeps calibration scripts,
tests, GLO-30 comparison reports, overlap replay checks and coverage evidence
separate. Both ZIPs contain checksums for their members. Extract both into the
same folder to run the optional calibration scripts; production needs only the
workflow archive plus the downloaded contour files.

```sh
python3 package_workflow.py --output-dir data/release/soi-dem-30m \
  --source-data-license ../../indianopenmaps/DATA_LICENSE.md
```

## Attribution and sources

Elevations: Survey of India source contour `VALUE` attributes from the
[soi-dem release](https://github.com/ramSeraph/indian_land_features/releases/tag/soi-dem).
Water polygons: Bhuvan/NRSC CARTODEM v3r1. Comparison and river-mask data:
[Copernicus DEM 2021 AWS distribution](https://copernicus-dem-30m.s3.amazonaws.com/readme.html).
WBM classes are defined in the
[Copernicus DEM handbook](https://dataspace.copernicus.eu/sites/default/files/media/files/2024-06/geo1988-copernicusdem-spe-002_producthandbook_i5.0.pdf).

The interpolation port and adapter are **GPL-3.0-or-later**, derived from
Christoph Hormann's Imagico implementation. `LICENSE.imagico.txt` and
`imagico_provenance.json` accompany the workflow and per-sheet audit packages.
The code license is not a claim to relicense upstream data. Release handling
uses [gh-release-tools](https://github.com/ramSeraph/release_tools).
The source release lists CC0 1.0 with requested attribution to Datameet and
Survey of India; its linked data-license text is copied into the workflow archive
as `DATA_LICENSE.source.md`. Retain that attribution with derived elevations.
