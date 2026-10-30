# National CartoDEM conversion and river handoff

Method record: 24 October 2026. This folder is also the local handoff for other
jobs. All paths below are relative to an arbitrary build workspace.

Current mask selection is **GLO-30 river class 3 throughout the DEM footprint**;
see [the source-change and execution record](GLO30_INTEGRATION.md). The WRIS
pilot and restart measurements below are historical and remain labelled as such.

## Agreed processing

1. Verify/download the pinned original CartoDEM sheets, EGM2008 grid and GLO-30
   classifications covering the full DEM footprint.
2. Mask original `-1024` and `-255` before computing `H = h - N`. Store both as
   Float32 `-255`, declared as nodata, on the original grid.
3. Remove GLO-30 ocean class 1. Every newly clipped pixel becomes the same
   `-255` nodata. No GLO-30 elevations are used.
4. Fit shared longitudinal profiles within GLO-30 river class 3. Lake class 2
   is excluded from fitting and bank samples. Heights are estimated from valid
   corrected CartoDEM bank samples in class 0. Preserve all source and ocean-mask
   nodata through this stage.
5. Export candidate native-grid Float32 COGs, then resample numeric heights for
   zooms 0–12 and encode quantized lossless WebP TerrainRGB and PMTiles.

Outside river class 3, retain CartoDEM's existing waterbody treatment apart from
datum conversion and ocean clipping. Stripe correction is not included.

All missing tile pixels use `(1, 134, 160, 254)`: a zero-metre placeholder with
alpha 254. Valid zero height uses alpha 255. The marker follows
[Iconem Terrain Viewer](https://terrain-viewer.iconem.com/docs/dev/demdiff-protocol/#holes-and-why-they-are-alpha-254).
A separate native-resolution validity mosaic prevents numeric COG overviews
from silently restoring holes. Pixel-centre coverage is sampled with nearest
neighbour; subpixel holes can disappear when no target centre samples them.

## Memory and sheet boundaries

Current geometry uses Whitebox raster cost allocation for water-to-centreline
ownership. The full GLO-30 context that previously exhausted memory completed in
25.83 seconds with 1,145.03 MiB peak observed process-tree RSS. The small skeleton
graphs remain in SciPy. See [ALLOCATION.md](ALLOCATION.md) for provenance, tie
differences, validation and cache replacement.

`staging.py` runs each full-coverage geometry extraction in a separate process.
The historical WRIS mode also isolates polygon rasterization. A long-lived WRIS
mask process was killed after 308 sheets at roughly
1.6 GiB resident anonymous memory. The next sheet succeeded by itself at about
402 MiB with a bounded GDAL cache. Process isolation releases retained GEOS/GDAL
allocations without changing mask membership or the geometry algorithm.

Geometry retains the existing 512-pixel halo. Each sheet owns its core samples
once; all skeleton nodes and water assignments use shared national pixel IDs.
Checks reject a halo owner absent from the assembled skeleton. Sheet borders
are not accepted as new outlets or independent fitting boundaries.

In WRIS comparison mode, the 45 invalid source polygons are handled only under the explicit native-grid
equivalence checks in [WRIS_INTEGRATION.md](WRIS_INTEGRATION.md). A source whose
normalization changes even one selected native pixel stops the build.

The current GLO-30 mask is already categorical raster data. It uses the verified
input masks directly; WRIS masks, geometry and sampled routing cannot be reused
after the source change. See [GLO30_INTEGRATION.md](GLO30_INTEGRATION.md).

## Direction evidence

`regional.py` runs WhiteboxTools 2.4.0 in overlapping contexts: one sheet plus
1024 native pixels on each side, clipped to the mosaic bounds, at most 40 million
cells. The stages are `BreachDepressions`, `D8Pointer`, then
`D8FlowAccumulation --pntr --out_type=cells`. The checksum-pinned precompiled
Linux ARM64 CLI uses one thread. A 20 ms RSS sampler terminates a tool observed
above 1100 MiB; this is not a kernel-enforced memory ceiling. Temporary rasters
are removed after a successful context and retained after failure for diagnosis.
Evidence, tool timings, observed memory use and raster hashes remain.

Only core-sheet river samples contribute, so overlap does not double-count
evidence. Exclude invalid directions, ocean, original holes and their uncertainty
buffer, the context-edge buffer, and nonpositive flow accumulation. Valid cells
beside NULL data or a context edge seed an uncertainty mask traced downstream
along the actual Whitebox pointers. Exclude these potentially incomplete
contributing areas; they do not establish an outlet. Store aggregated direction
vectors, weights and sample counts per
global skeleton owner. Small-hole interpolation is confined to the disposable
routing surface; it never fills published elevation rasters. Breaching also
affects only this surface. Whitebox's TIFF input uses predictor 1; every sample,
the grid and nodata are verified identical to the original predictor-3 analysis
raster before routing. Output grids and NULL footprints are checked again.

This is supporting direction evidence, not a complete national catchment
delineation. Large rivers entering a context may lose routing evidence. The
existing bank-trend and joint-profile logic then decides support. Independently
supported directions remain distinguishable from joint-profile hypotheses.

The first national routing context failed in GRASS after geometry for all 543
sheets finished. The independent [Whitebox routing and watershed trial](WHITEBOX_TRIAL.md)
passed: precompiled and source-built tools produced identical rasters. Whitebox
is now the national runner's default. GRASS is only an optional historical
comparison backend and is not needed by this path. Its intermittent failure's
cause remains unresolved. `Watershed` validates a trial catchment against
accumulation; production needs only the three routing stages above.

## Connected fitting and projection

`national.py` builds a sparse graph of all skeleton nodes and adds connectivity
from valid native water contacts, including contacts across sheet borders.
Independent connected components are staged in `components/chunks.sqlite` and
fitted one at a time with the existing `joint.py` solver. This reduces memory
without fitting adjoining portions independently. The largest connected
component still determines peak fitting memory.

Each junction has one fitted elevation shared by its adjoining reaches. The
existing quantile, station spacing, bank-support, direction and adjustment
checks remain unchanged. Raw water-pixel heights are not exact boundary pins.
Entire model groups remain unchanged when supported and retained surfaces would
otherwise meet without a shared fit. Unsupported components retain their
corrected CartoDEM elevations. This is the current complete-component rule,
not the rejected per-pixel source-pin experiment.

`components/models.sqlite` checkpoints completed fits, component decisions,
station values and native skeleton projections. Recipe hashes prevent reuse
after source/code/settings changes. Water pixels receive their existing
in-channel nearest-skeleton projection; no shortcut crosses land or islands.

Application checks exact native COG values and validity, unchanged nonriver and
nodata pixels, ocean nodata, downstream station order, and treatment interfaces
within and across sheets. Per-sheet audit reports record introduced native
adjacent-water gradients and their largest examples.

## Closed-loop fitting restart, 25 October 2026

Geometry and regional routing completed for all 543 sheets. The first national
fit stopped after 6,358 of 52,807 connected components, at component 6358 in
`cdnh43m.tif`. Its four skeleton pixels form a closed loop with no junction or
terminal station. An empty station-index array defaulted to floating point,
causing a NumPy `IndexError` before the existing unsupported-component handling.

`compress_network` now explicitly creates integer station indices, including
the empty case. Closed loops still have no fitted surface and retain their
corrected source heights. The fitting parameters and direction rules are
unchanged. Both the exact failing four-node geometry and a two-sheet COG
integration test reproduced the failure before the fix; all 62 tests passed
afterwards, including unchanged loop heights and nodata in exported COGs.

The change affects fitting only; geometry extraction and routing implementations
are unchanged. Preserve the old partial `components/models.sqlite`, its
`fit-status.json`, and the failed build status before restarting. Recompute the
partial fit under the new source hash, while retaining the verified partition,
geometry and routing checkpoints and their original recipes. Do not rewrite
checkpoint hashes to accept changed code. The portable
[restart record](pilots/fit-empty-loop-v1/restart.json) records the preserved
files, source snapshot and validation. The restart command is:

```bash
python run_build.py --workdir /path/to/build --river-source glo30 \
  --resume-from fit --detach
```

This records a fitting restart, not completed national COGs or tiles.

## Quality status

Outputs remain candidates until the recorded gradients are reviewed. The WRIS
joint-fit pilot previously showed up to about 9.59 m of introduced adjacent-pixel
change within a fitted profile, despite shared junctions and no partial-treatment
interfaces. Short model segments and projection within wide channels are known
limitations. The runner does not claim those have been fixed. The historical
10 m statistic is diagnostic, not an acceptable-error standard or a new gate.
Candidate files are built locally in the workspace; the runner does not publish
them to an external release.

## Historical WRIS validation and kickoff record

All 543 corrected coastal sheets were verified before the national restart.
The ocean mask removed 26,380,306 previously valid samples; all became `-255`.
The source contained 117,449,248 samples of `-1024` and 375,500,287 of `-255`.
All original holes were retained, including 61,371,759 on GLO-30 land. There are
6,517,950,159 valid prepared samples after ocean masking.

The 51-test build suite passed, including native mask precedence, both nodata
codes, alpha/quantization round trips, COG preservation, cross-sheet profiles,
component equivalence and isolated-worker geometry. The four tile integration
tests also passed after adding candidate status and WRIS attribution to archives.

The real `pilot-national-v1` completed geometry, two regional GRASS contexts,
4,399 staged components and both COG exports. It treated 106,545 of 481,308 valid
river pixels (22.14%) and preserved 646,658 nodata samples. The native shared-fit
audit reported no partial-treatment interfaces or downstream station violations.
The maximum introduced adjacent-water jump remained 9.5933837890625 m; that
quality issue is not resolved by this execution work.

Changing the GRASS contexts changed 52 fitted water pixels in D43V by at most
3.015594482421875 m; D43W was identical. Validity was identical in both sheets.
To separate this from batching effects, `validate_partitioning.py` refitted the
whole two-sheet network with exactly the new routing evidence. All 38,960 treated
skeleton nodes had identical treatment and **exactly equal Float32 heights**
to the component-by-component fit. The geometry arrays also matched the earlier
pilot exactly. The fitting method was preserved; context-dependent direction
evidence accounts for the changed pilot elevations.

The portable pilot records are in [pilots/national-v1](pilots/national-v1/),
including [the equivalence result](pilots/national-v1/partition-equivalence.json).
The initial national restart used `run_build.py --resume-from geometry --detach`,
reusing verified prepared elevations and masks. After the GRASS failure, the
Whitebox restart at 08:42:22 UTC on 24 October 2026 used
`run_build.py --resume-from routing --detach`. Its source snapshot is
`reports/source-20261024T084222Z-whitebox.tar.gz`, SHA-256
`90723fa152391f8c94da801f53fd6e22f6b6d34063aa6a89dfdb4561185ac875`.
The previous failed status is preserved separately.

The 55-test suite passed with the Whitebox adapter. The first production context
completed in 43.64 seconds and exactly matched the independent precompiled trial:
122,860 accepted samples and all direction vectors, weights and counts for
31,358 owners. The accumulation TIFF also has the identical SHA-256. See
[the production equivalence record](pilots/whitebox-routing/production-equivalence.json)
and [restart record](pilots/whitebox-routing/restart.json).

The first four contexts completed in 43.64, 79.16, 78.67 and 24.35 seconds:
56.45 seconds per sheet on average. Extrapolating this small coastal sample gives
about 8.5 hours for 543 routing contexts, not for the complete build; inland
coverage and complexity can change that rate. Fitting, COG export, TerrainRGB
and PMTiles follow. `reports/build-final.json` and `logs/build-final.log` contain
live progress. This kickoff record does not claim national completion.

## Reproduction and restart

From the `egm2008` source directory, install the pinned Python dependencies
using Python 3.14. On Linux ARM64, install the precompiled tools and run:

```bash
.venv/bin/python install_whitebox_binary.py --workdir /path/to/build
.venv/bin/python run_build.py --workdir /path/to/build --detach
```

This includes downloads. Existing original files are verified before reuse.
`inputs.lock.json` and `wris.lock.json` pin input checksums. No execution-host
connection information is required by the scripts.

For a bounded validation of the same regional/component path:

```bash
.venv/bin/python run_build.py --workdir /path/to/build \
  --label pilot-whitebox-v1 --sheets cdnd43v.tif cdnd43w.tif \
  --resume-from geometry --through cogs --detach
```

This assumes the two coastal sheets already passed `build.py prepare`. The
small-grid path in `build.py` remains available for historical comparisons.

The coordinator writes `logs/build-LABEL.log` and `reports/build-LABEL.json`.
Resume from a failed stage with `--resume-from geometry` (or `routing`, `fit`,
`cogs`, `mosaic`, `tiles`, `package`). It stops on errors; it never silently skips
a failed component or audit. Use a fresh label when a recipe changes after
geometry/fitting has been cached. Completed per-sheet and per-component work is
verified and reused under an identical recipe.

Key records under `rivers/LABEL/` are `geometry-recipe.json`, `geometry.json`,
`routing/report.json`, `components/partition.json`, `components/models.sqlite`,
`fit.json`, `cogs.json`, and `audit.json`. The COGs are in `output/LABEL/cogs/`.
WRIS comparison mode additionally writes `wris-masks.json`.
The full `final` build continues to tile databases and PMTiles; after packaging,
`reports/RUN_NOTES.md` records counts, hashes, versions and the candidate status.
