# Experimental river flattening

This is the method used to produce the final CartoDEM COGs. The input elevations
have already undergone EGM2008 correction and GLO-30 ocean masking. The method
uses **GLO-30 WBM river class 3 for geometry and CartoDEM for heights**. It does
not copy GLO-30 elevations, flatten lakes, or fill nodata in the published DEM.
See [the reproduction commands](../README.md) for the complete pipeline.

## Geometry and bank observations

`production.py` verifies the prepared elevation sheets and their WBM hashes,
then builds shared native-grid elevation and classification mosaics. WBM class
0 is land, 1 is ocean, 2 is lake and 3 is river. Class 255 denotes missing
classification coverage in the working mosaic and cannot support a valid
CartoDEM pixel.

`geometry.py` works directly on the class-3 raster; no river polygon extraction
or GDAL polygonization is required. Each sheet is read with a 512-pixel halo.
The river mask is skeletonized, short bank-corner branches are pruned, and
degree-two chains are compressed into reaches. The pruning threshold is the
larger of 90 m and twice the junction's distance from the bank. Junction pixels
share a station, and retained reaches receive stations at a target spacing of
180 m. Pruning changes the centreline representation, not the river footprint.

Every water pixel is assigned to a retained skeleton pixel through the river
mask. `allocation.py` uses WhiteboxTools **CostDistance** followed by
**CostAllocation**, with unit cost and eight neighbours. Land and islands are
nodata barriers. A disposable metric grid uses the window's constant east-west
and north-south pixel distances; published grids are never reprojected through
that grid. This avoids the memory cost of explicitly constructing a graph with
one node for every water pixel. Seed identity and water-mask coverage are checked
after allocation. Equal-distance ties follow the pinned Whitebox implementation.

Bank observations are valid corrected CartoDEM samples in GLO-30 class 0 within
two native pixels of the river. Lake, ocean, river and missing samples are
excluded. Bank samples enter the network through their nearest river pixel and
its in-river skeleton assignment. No synthetic water depth is subtracted.

Each core sheet contributes its samples once. Nodes use global native-grid
pixel IDs. Native water contacts, including contacts across sheet boundaries,
connect the assembled skeleton graph. A sheet boundary is not an outlet or an
independent fitting boundary. Geometry extraction runs in isolated processes
for the full inventory so retained GDAL and array allocations are released.

## Downstream direction evidence

`regional.py` runs overlapping WhiteboxTools 2.4.0 contexts: one sheet plus a
1024-pixel halo on each side, limited to 40 million cells. The commands are:

```text
BreachDepressions
D8Pointer
D8FlowAccumulation --pntr --out_type=cells
```

These operate on a **disposable routing surface**. Small internal holes may be
interpolated within a 64-pixel search distance for routing, and depressions are
breached there. Neither interpolation nor breaching changes the elevation
baseline or final COGs. Whitebox's TIFF input uses predictor 1; its samples,
nodata and grid are checked against the analysis raster before use.

Routing evidence excludes original invalid samples, an eight-pixel uncertainty
buffer, the context-edge buffer, and invalid/nonpositive accumulations. Possible
missing upstream area at context boundaries and nodata edges is traced
downstream using the actual D8 pointers and excluded. Only core-sheet river
samples contribute to the evidence. Direction vectors are accumulated per
global skeleton owner, weighted by `min(log(1 + accumulation), 12)`.

`directions.py` combines that evidence with:

- Robust CartoDEM bank trends, estimated by Theil–Sen slope. A bank trend needs
  at least three observations over 90 m, a 90% slope interval excluding zero,
  and at least 3 m of inferred drop.
- Ocean-connected terminal mouths, allowing more than one mouth in a network.
- Outlet-path suggestions from the selected RivGraph functions.
- Continuity at junctions, where accepted directions support the remaining link.

Terrain vectors are projected onto the centreline tangent. A reach needs at
least three usable station votes and 75% agreement. A routing direction without
bank evidence also needs outlet-path agreement. Conflicting bank/routing
evidence, unexplained interior sources/sinks and directed cycles are left
unresolved at this stage. A shortest-path suggestion alone does not establish
flow direction.

This provides local direction evidence, not a national watershed delineation.
The production path does not run Whitebox `Watershed`; that tool was evaluated
in an earlier trial. GRASS is not required by the production commands.

## A joint profile across connected reaches

`national.py` assembles the complete skeleton graph, finds its connected
components, and stages them in SQLite. Each component is fitted once with
`joint.py`, then projected into all of its sheets. Components are independent
work units; adjoining portions of a river are not fitted independently merely
because they lie in different TIFFs.

All reaches meeting at a junction share **one elevation unknown**. Raw river
pixels are not exact elevation pins. `profiles.py` uses SciPy sparse matrices
and Clarabel to minimize weighted quantile residuals against bank observations
plus a squared-curvature penalty along the longitudinal profile:

```text
minimize  sum(weight * quantile_loss(bank_height - station_height))
          + curvature_penalty

subject to station heights not increasing downstream
           and supported station heights <= their bank upper bound
```

Distances in the curvature calculation are in metres and may be unequal.
Stations receive equal total observation weight, so a station with more bank
pixels does not dominate solely through its sample count. The recorded settings
are:

| Setting | Value used |
| --- | --- |
| Station-spacing target | 180 m |
| Quantile and supported-station bank upper bound | 20th percentile |
| Curvature stiffness | 2.0 |
| Bank observations to support a station | At least 4 |
| Support within a connected model group | At least 3 stations and at least 60% of stations |
| Maximum along-network distance to support | 1,200 m |
| Allowed fitted departure from supported bank quantiles | 95th percentile ≤15 m; maximum ≤40 m |

These are experimental parameter choices, not DEM accuracy guarantees or a
general hydro-flattening standard.

The first fit constrains independently supported directions. Where direction
is unresolved, its endpoint heights propose a direction for a second fit.
Endpoints within 1 mm propose a level reach. Equal-level relations and strongly
connected direction groups share height unknowns, and the second fit enforces
the resulting non-increasing or level profiles. These additional directions
are recorded as **joint profile hypotheses**, distinct from independent evidence.
Every expanded reach is checked again after fitting and Float32 conversion.

Groups with insufficient bank support, excessive gaps or adjustments, unusable
geometry, or failed numerical checks keep the corrected source heights. Closed
loops without a modelled reach are also retained. Negative elevations are
allowed; a valid value equal to the reserved `-255` sentinel fails.

`complete_components()` checks native water contacts. If a fitted group would
touch an unresolved group or a separately fitted group without a shared model,
the affected fitted group is retained unchanged. This is repeated until the
partial-treatment interfaces disappear. The method does not accept a 10 m
boundary discontinuity or pin arbitrary source water pixels to make a fit pass.

## Apply the surface and export COGs

Fitted station heights are interpolated along the skeleton. Each river pixel
receives its owner's interpolated height only if it has a supported model and
was valid in the corrected ocean-clipped baseline. All original nodata and all
ocean-clipped pixels remain `-255`. Unsupported river pixels and every nonriver
pixel retain their baseline values exactly.

Pixels assigned to the same skeleton position share a height. In bends or wide
channels this discrete allocation is not a surveyed cross-section; it can
produce raster gradients even when the longitudinal station profile is monotone.
COG export and the independent native audit check pixel values, validity,
nonriver preservation, downstream station order, treatment interfaces, and
sheet-border contacts. Averaged COG overviews are built after native values have
been finalized; the tiler uses a separate native validity mosaic.

## Results and limitations of this run

The full build produced 543 COGs, with:

- **6,388,041 of 42,457,442 valid river-mask pixels treated**, about 15.05%.
  Remaining river pixels retain their corrected, ocean-clipped baseline.
- 52,807 connected skeleton components assembled for fitting.
- No recorded downstream violations in the stored fitted station profiles.
- No partial-treatment native-water interfaces, including across sheet borders.
- A **maximum 33.9862 m increase in an adjacent-water height step** relative to
  the baseline. This is `abs(after_b-after_a) - abs(before_b-before_a)` for
  neighbouring valid river pixels; it is not a measured elevation error against
  surveyed truth.

The last result is why shared junctions and monotone station profiles must not
be described as proof that every resulting river surface is smooth or correct.
Native projection and short model segments can still create steps. No stripe
correction, surveyed river-height calibration, complete watershed solution or
hydrological validation is claimed. The source manifests retain the build's
candidate/gradient-review status, and the release describes the treatment as
experimental.

Runtime records are `rivers/final/geometry.json`, `routing/report.json`,
`components/partition.json`, `components/models.sqlite`, `fit.json`, `cogs.json`
and `audit.json`. Recipes contain input, tool and source-code hashes. Temporary
allocation/routing rasters are deleted after successful use; geometry,
aggregated direction evidence, fitted models and reports remain available for
resuming later stages.

## Method references and reused code

[Schwanghart and Scherler (2017)](https://doi.org/10.5194/esurf-5-821-2017)
describe quantile-based longitudinal profile correction and smoothing. The
solver here is an independent implementation adapted to bank observations and
shared river junctions; it does not execute TopoToolbox or reproduce the
WorldDEM editing workflow.

The limited RivGraph code used for direction suggestions is vendored from
[mwort/RivGraph at revision 0a0a713](https://github.com/mwort/RivGraph/tree/0a0a7135752894f451344bbf3521cc0f85f6f426).
[vendor/provenance.json](vendor/provenance.json) lists the selected functions,
source hashes and local edits; [the BSD-3-Clause license](vendor/LICENSE-RivGraph.txt)
is included. WhiteboxTools is downloaded as a pinned standalone binary by
`install_whitebox_binary.py`; its package licenses and binary hash are retained
with the installation.

Earlier GLO-30, WRIS, GRASS and source-pin experiments are archived under
[old/egm2008/river_flattening](../../old/egm2008/river_flattening/README.md).
Their treatment percentages and transition statistics do not describe this
national GLO-30/Whitebox run.
