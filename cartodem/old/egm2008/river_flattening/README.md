# River flattening method and implementation record

This folder is the local handoff for other jobs that need the CartoDEM river
treatment. It is part of the reproducible build source, independent of where the
large raster processing runs.

**Current method: GLO-30 WBM river class 3 and shared connected profiles.**
[NATIONAL_BUILD.md](NATIONAL_BUILD.md) describes the resumable national runner,
memory isolation, global component assembly and candidate-output review.
See [GLO30_INTEGRATION.md](GLO30_INTEGRATION.md) for coverage beyond India,
mask selection and the restart record. [WRIS_INTEGRATION.md](WRIS_INTEGRATION.md)
records the earlier polygon comparison and development of the shared-junction
fit. The following GLO-30 figures describe historical direction pilots, not the
current GLO-30/Whitebox run.

Water-to-centreline assignment now uses Whitebox `CostDistance` and
`CostAllocation`, avoiding the large water-pixel graph that exhausted memory.
[ALLOCATION.md](ALLOCATION.md) records the implementation, complete-worker memory
test and measured differences against the preceding two-sheet pilot.

`profiles.py`
implements the shared longitudinal profile solver; `geometry.py` constructs a
centreline network and assigns water pixels; `study.py` runs a bounded real-data
pilot. `production.py` now assembles shared networks across sheets and exports
COGs. The GRASS/RivGraph integration passes 19 river tests (33 tests across the
complete build). Its latest full D43V/D43W pilot treats 61,291 of 160,151 valid
river pixels (38.3%), after retaining ten groups that introduced large transitions
against unchanged river elevations. Independent pixel, validity and downstream
checks pass. See [the integration method](DIRECTION_INTEGRATION.md).
The bounded D43V v4 pilot treated 15,108 of 17,381 valid
river pixels, but its crop hid additional ocean outlets. The wider D43V/D43W test
treated only 27,218 of 160,151 valid river pixels (17.0%) under the current
earlier conservative rules. Land, lake class 2 and original nodata remain unchanged. See
[the pilot findings](PILOT_RESULTS.md) and [its complete report](pilots/d43v-v4/report.json).
Two updated pilot COGs have been produced and checked; no nationwide river product has
been produced or published. The bounded pilot's 86.9% coverage must not be used
as an estimate of nationwide coverage.

[Existing-code reuse review](REUSE_REVIEW.md) records GRASS, RivGraph and GRIT
implementations, exact revisions, licenses and integration gaps. `production.py`
uses the scoped adapter in `directions.py`. The current national runner obtains
Whitebox evidence through `regional.py`; `routing.py` retains the historical
small-grid GRASS path.
The old `study.py` retains the original harmonic method for historical comparison.

## Scope and inputs

- Elevations come exclusively from CartoDEM, after EGM2008 correction.
- GLO-30 WBM class **1** removes ocean coverage and class **3** supplies river
  geometry throughout the DEM footprint. No country boundary clips the mask.
- Lake class **2** is excluded from river treatment and bank observations.
  There is no new lake-flattening stage. Elevations outside class 3 retain their existing
  CartoDEM treatment, apart from the datum correction and ocean masking.
- Neither GLO-30 elevations nor invented river depths are inputs.
- Original `-1024` and `-255` are masked before datum conversion. River processing
  preserves those holes and every newly clipped GLO-30 ocean pixel. Final COGs
  store all of them as `-255`, with `NoData=-255`.
- Stripe correction is outside this build's scope.

## Basis for the method

[USGS DEM surface-treatment guidance](https://www.usgs.gov/ngp-standards-and-specifications/lidar-base-specification-digital-elevation-model-surface)
describes river surfaces that are level across the channel and descend along it.
It also requires consistency with surrounding terrain. This is guidance for the
desired surface; it does not supply a complete automatic algorithm for noisy
30 m CartoDEM measurements.

The [GRASS hydro-flattening workflow](https://grass-tutorials.osgeo.org/content/tutorials/hydro_flattening/hydro_flatenning.html#adding-break-lines)
uses river subdivisions and bank samples. Its lesson here is to estimate water
levels locally along a river, rather than assigning one height to an entire
connected river system. Its lidar-specific thresholds are not adopted unchanged.

[Schwanghart and Scherler (2017)](https://doi.org/10.5194/esurf-5-821-2017)
describe constrained quantile smoothing for longitudinal river profiles. The
[TopoToolbox implementation documentation](https://topotoolbox.github.io/pytopotoolbox/_modules/topotoolbox/stream_object.html#StreamObject.crs)
provides a working reference for that formulation. Our numerical kernel is an
independent implementation using NumPy, SciPy and Clarabel; it does not run or vendor
TopoToolbox.

The [NGA pydrodem implementation](https://github.com/ngageoint/pydrodem/blob/main/models/burn_rivers.py)
was also examined. Its complete hydrological-conditioning workflow includes
stream burning and other changes beyond this cartographic task, so that entire
workflow is not being applied.

### How WorldDEM / GLO-30 handles rivers

The [Airbus production paper, section 3](https://isprs-archives.copernicus.org/articles/XL-7-W3/1149/2015/)
describes automated water detection from radar amplitude imagery, operator
correction and classification, and shoreline-derived river elevations descending
in 0.5 m steps. Operators participate in selecting river heights. It does not
give a complete automatic algorithm for resolving multiple outlets. Our
single-outlet restriction is a prototype policy, not the documented WorldDEM
method.

GLO-30 inherits edited WorldDEM data. The
[Copernicus handbook, page 15](https://dataspace.copernicus.eu/sites/default/files/media/files/2024-06/geo1988-copernicusdem-spe-002_producthandbook_i5.0.pdf)
warns that resampling can aggregate nearby hydrographic features and that
hydrological consistency is not guaranteed. Published WBM connectivity alone
therefore cannot be assumed to reproduce the original edited network.

## Shared profile solver: implemented

Represent the river network with stations and directed edges, ordered upstream
to downstream. All reaches meeting at a junction use **the same station index**.
Distances along reaches are in metres; observations are corrected CartoDEM
elevations in metres. Sheet identifiers are irrelevant to the fit.

For an observation `z` assigned to station `i`, fit a water elevation `h[i]` by
minimizing asymmetric quantile residuals plus squared profile curvature:

```text
sum(weight * quantile_loss(z - h[i])) + sum((s * second_derivative(h)) ** 2)

subject to h[upstream] >= h[downstream] on every edge
and optional h[i] <= bank_upper_bound[i]
```

Second derivatives use the actual, possibly unequal station spacing. Junctions
share elevations; separate reaches may have different slopes. Split and rejoined
channels can share upstream and downstream stations, provided their directions
form a directed acyclic graph.

Each sampled station receives equal total observation weight by default, so
stations with more bank pixels do not dominate solely because of pixel count.
The current **trial** parameters are quantile `0.2` and stiffness `2.0`; these are
not yet selected production values. The smoothing scale follows station spacing
and network size, so the dimensionless stiffness alone is not a distance in
metres. The current pilot uses a 180 m station-spacing target, land samples within
two native pixels of the river and the per-station 20th bank percentile as an
upper bound. These values still require broader validation.

The solver rejects cycles, inconsistent reach directions, nonfinite observations,
unobserved components and failed convergence. Solver inequality tolerance must
be below 1 mm; an adjustment of at most 1 mm removes numerical uphill steps before
Float32 storage. Larger violations fail. Negative elevations are allowed. A valid
height exactly equal to the reserved output sentinel `-255` fails explicitly.

## Geometry, flow direction and cross-sheet work

A water mask supplies a footprint, not downstream direction. Estuaries, braided
channels, disconnected mask fragments and noisy bank heights can make direction
ambiguous. `fit_network` deliberately **requires** directions; it does not silently
infer them from one low pixel. The pilot's direction candidates and rejection
reasons are explicit in its report. Unsupported components remain unchanged.

`geometry.py` skeletonizes the selected river class, removes short branches toward bank corners,
and compresses degree-two chains into stations. Adjacent junction pixels share
one station. Pruning changes the centreline representation, never the water
footprint. The pruning threshold is the larger of 90 m and twice the junction's
distance to the river edge; window-edge ports and standalone channels are kept.

Each water pixel is assigned to its nearest retained skeleton pixel by shortest
distance **inside the river mask**. This avoids a shortcut across a meander or
island. Bank samples enter through their nearest river pixel, then use that
pixel's skeleton assignment. Valid corrected CartoDEM samples outside the river
and ocean supply bank observations. Current GLO-30 mode uses class 0 only,
excluding lake class 2 as well as nodata from bank observations.

The historical bounded pilot chooses an outlet only when one terminal adjoins ocean class 1, or when
the lowest terminal bank quantile is at least 3 m below the next candidate. The
latter is a hypothesis requiring validation, not an observed flow measurement.
Multiple ocean outlets were unsupported by that historical method. The integrated
driver permits multiple mouths and evaluates each reach; see
[direction decisions](DIRECTION_INTEGRATION.md#direction-decisions).

Given an outlet, a harmonic potential with outlet=0 and terminal inlets=1 supplies
a geometric direction candidate. These are **dimensionless direction values**,
not water heights. Unlike using distance to the outlet, this permits both arms of
an unequal-length braid to point from their common upstream junction to their
common downstream junction. Tied potentials, closed loops without usable ends,
and inconsistent networks remain unsupported. This geometric routing is our
prototype choice, not a claim that the cited CRS paper establishes direction.

Bank-height contradictions over 5 m reject a direction inferred solely from
bank levels. With a single coast-connected outlet, contradictions are retained
as diagnostics: noisy bank measurements are one of the reasons for fitting the
surface in the first place. The fit must still pass support and adjustment
checks. At least 60% of stations need four bank samples; the bounded study defaults
to a maximum gap of 600 m from supported stations. The v4 pilot explicitly tested 1,200 m to span
missing bank observations while retaining every original raster hole. Fitted
station adjustments must have a 95th percentile at most 15 m and a maximum at
most 40 m relative to supported bank quantiles. The production driver currently
defaults to the same 1,200 m gap tested in v4 and the two-sheet run; that setting
has not been validated nationwide. These are trial thresholds.

Neighbouring sheets are needed for bank observations and channel geometry.
A halo supports local sampling but does not establish a shared profile for a
long river. Assemble connected reaches across sheet boundaries, fit each network
once, and then project the same station heights back into each native sheet.
Independent fits per sheet are not acceptable. The geometry API already returns
common-grid pixel IDs, and a test verifies identical ownership across two
overlapping windows. `production.py` implements per-sheet geometry caches, shared
network assembly and fitting, and native-grid COG application. A two-sheet real
pilot exercised those stages. Assembly and fitting across all 543 source sheets
have not yet run; a bounded pilot is not a substitute for that validation.

## Raster application and validation

The pilot applies the interpolated profile only to pixels with valid CartoDEM elevations
inside GLO-30 river class 3. Keep nonriver elevations, lake class 2 and original
nodata unchanged.
Pixels assigned to the same skeleton position receive the same interpolated
height. At bends this is a discrete geodesic assignment, not a surveyed transect;
cross-section validation is still needed beyond the synthetic straight-channel
test. Arbitrary two-dimensional smoothing is not a substitute for the profile.

Further pilots must check geometry, direction, station support, bank-to-bank levelness,
downstream monotonicity, confluences and sheet seams. Inspect changes against the
coastal-only EGM2008 baseline and report unsupported reaches and large residuals.
After those checks, write publishable native-grid Float32 COGs with numeric
overviews and declared `NoData=-255`. TerrainRGB generation follows those final
elevations, with the existing quantization and alpha-254 nodata convention.

## Reuse and tests

The numerical entry point is `river_flattening.profiles.fit_network`. It returns
one Float32 elevation per shared station and a report of numerical checks.
It does not read rasters, manufacture directions, fill nodata or flatten lakes.

From the parent `egm2008` directory:

```bash
python -m unittest -v river_flattening.test_profiles river_flattening.test_geometry river_flattening.test_production river_flattening.test_directions river_flattening.test_boundaries
```

Tests exercise confluences, a split-and-rejoined channel, an isolated low outlier,
missing bank support, negative heights, bank caps, reversed reaches and cycles.
The v4 run passed the nine numerical/geometry tests; the later production test
adds a synthetic two-sheet integration case. The v4 report also
records hashes of its input windows and a snapshot of the processing code under
`pilots/d43v-v4/code/`. Earlier v1–v3 directories are exploratory diagnostics;
v4 is the first pilot with a complete source snapshot and input-window hashes.
