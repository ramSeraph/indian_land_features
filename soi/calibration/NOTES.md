# Historical calibration notes

Run the commands below from `soi/`. These notes retain earlier experiments; some commands and outputs have since been replaced. Use [the current workflow](../PRODUCTION.md) for production.

# SOI contour-to-DEM pilot

For the reproducible local batch, pinned downloads, production/calibration
separation, release uploads and current validation results, see
[PRODUCTION.md](../PRODUCTION.md). Independent-sheet outputs retain their original
constraints for later neighbour reconciliation and water conditioning.

This implements a single-sheet baseline for the first part of
[issue #2](https://github.com/ramSeraph/indian_land_features/issues/2), using the
released SOI contours and their `VALUE` elevations. It does not use the coloured
WMS DEM images or the 1:25,000 NHP maps. The earlier WMS experiments remain in
`download_wms_images.py` and `rgb_to_dem.py`.

The default engine is a Python adaptation of
[Imagico's contour interpolation method](https://www.imagico.de/pov/earth_interpolate.php),
using NumPy and SciPy without CImg. It runs directional interpolation,
unknown-height assignment, curvature refinement, extremum constraints,
gradient paths, and final interpolation/refinement. The previous SciPy linear
TIN is available with `--interpolation linear` for historical comparisons.

The numerical port in `imagico_interpolate.py` and raster adapter in
`imagico_dem.py` use **GPL-3.0-or-later**, matching Christoph Hormann's original
code. [LICENSE.imagico.txt](../LICENSE.imagico.txt) contains the full license;
[imagico_provenance.json](../imagico_provenance.json) records the author, upstream
version, source URLs and hashes, and intentional changes. Original source is
retained in `contour_interpolate-0.3/`. SciPy filtering, explicit missing-data
masks and constraint protection mean this is not a bit-for-bit CImg reproduction.
Engine reports include hashes of the port, adapter, provenance and license.

## Local contours and sheet inventory

All four released `SOI_Contours_raw.*.parquet` partitions and their metadata
were downloaded on 2026-10-23 to `data/contours/soi-raw/`. The Parquets occupy
4.518 GiB and contain 7,362,891 source rows. `download_manifest.json` in that
directory records release URLs, sizes, SHA256 digests and successful verification.

`inventory_contour_sheets.py` streams those files and intersects their actual
geometries with all 5,113 WGS84 polygons in `index_50k.geojson`. It distinguishes
positive-length linework beyond a 2 m inward edge band from edge-only linework
and point contacts. It does not infer sheet coverage from feature bounding boxes.

```sh
uv run inventory_contour_sheets.py \
  data/contours/soi-raw/SOI_Contours_raw.parquet.meta.json \
  --output-dir data/contours/soi-raw/sheet-inventory --plot
```

The output directory must be new or empty. The completed inventory is in
`data/contours/soi-raw/sheet-inventory/`:

| Status | Sheets | Sheet list |
| --- | ---: | --- |
| Interior linework beyond the 2 m edge band | 2,321 | `sheets_interior.txt` |
| Linework confined to the edge band or boundary | 83 | `sheets_edge_only.txt` |
| Point contacts only | 21 | `sheets_touch_only.txt` |
| No intersecting usable linework | 2,688 | `sheets_no_data.txt` |

`sheet_coverage.csv` records the old and OSM sheet IDs, status and source-row
counts. `sheet_coverage.geojson` carries the same properties on sheet polygons;
`sheet_coverage.png` shows the coverage map. `report.json` records inputs,
definitions and scan statistics. The scan took about 49 seconds on this machine.

These are raw-data presence counts, including circular and triangular height
symbols. Rows are not deduplicated and can intersect multiple sheets. The 3,631
zero-length geometries are counted separately and excluded from coverage; all
7,359,260 remaining source rows intersect the index. Presence does not establish
complete coverage or usable elevations: 31 of the interior sheets have only
1–10 interior source rows. The inventory neither changes nor filters the DEM
pipeline's source constraints.

## Extract a sheet with iomaps

Use [iomaps](https://github.com/ramSeraph/indianopenmaps/tree/main/python) from the
local checkout containing the bounding-box pushdown and nested-attribute export
fixes. These fixes have not been published to PyPI as part of this work.

To extract the contours and water mask and build the DEM in one command:

```sh
uv run extract_sheet_dem.py --sheet 56M_9 \
  --iomaps ../../indianopenmaps/python/.venv/bin/iomaps \
  --output-dir data/contour_dem/56M_9
```

This command needs no toposheet, manual label reading, or review file. It uses
`index_50k.geojson` (downloading the default index if absent), caches the iomaps
extracts in `data/contours/<sheet>`, and records extraction provenance in
`extraction.json`. Use `--reuse-extracts` for existing local extracts; files with
cache manifests must match their recorded sheet, source, geometry and checksum.
Older extracts without manifests are explicitly recorded as having no extraction
manifest. The checkout now uses DuckDB 2.0's native asynchronous HTTP range
reader with connection pooling. In live comparisons, scan/export took 22 seconds
for 45F_4 and 60 seconds for 73D_15, versus roughly 12 minutes each with the
earlier DuckDB 1.5.6 setup. The resulting features, geometries and source values
match exactly. These are individual runs, not controlled benchmarks; metadata
requests add about 35 seconds to each 2.0 run.

`--http-client curl` or `--http-client httplib` can override the native selection.
The tested default needs no community HTTP extension. `cache_httpfs` and
`curl_httpfs` were not yet available for DuckDB 2.0 on macOS ARM64 when checked
on 2026-10-22. The complete Parquets are now available locally as described above;
this extraction command still uses the release URLs. Empty contour
extracts and incomplete downloads are rejected. An empty water extract is valid.
`--no-water` omits the optional smoothing step.

Toposheet GeoTIFFs are used only by the separate calibration workflow below.

For extraction and DEM construction as separate steps:

From this directory:

```sh
mkdir -p data/contours/56M_5
../../indianopenmaps/python/.venv/bin/iomaps cli extract \
  -s /dem/contours/raw/soi/ \
  -f index_50k.geojson --pick-filter-feature-kv id=56M_5 \
  -o data/contours/56M_5/contours.geojsonl

../../indianopenmaps/python/.venv/bin/iomaps cli extract \
  -s /dem/water-bodies/cartodem-v3r1/bhuvan/ \
  -f index_50k.geojson --pick-filter-feature-kv id=56M_5 \
  -o data/contours/56M_5/water.fixed.geojsonl
```

`iomaps` reads release metadata, selects every intersecting Parquet partition,
and filters features by the selected sheet polygon. For 56M_5 this selects
`SOI_Contours_raw.000.parquet` and `.010.parquet`. It retains complete intersecting
geometries, allowing circular symbols to be identified before clipping.

The published partition extents do not cover 47H_8 (16–16.25°N); their southern
limit is 18°N. An empty extract is not a zero-elevation surface.

## Build the baseline

```sh
uv run contours_to_dem.py data/contours/56M_5/contours.geojsonl \
  --sheet 56M_5 --resolution 30 \
  --water data/contours/56M_5/water.fixed.geojsonl \
  --output-dir data/contour_dem/56M_5-auto-symbols
```

The script chooses the sheet's UTM zone, identifies circular and triangular
symbols, and converts them into centre points. The circle test requires
a simple closed ring, radius at most 13 m, radial coefficient of variation at
most 3%, and roundness `4π × area / perimeter²` at least 0.98. Roundness rejects
arcs closed with a chord even when all vertices lie exactly on a circle.
`--spot-radius 0` disables automatic circle conversion.

Triangle detection assumes the same symbol scale as the GeoTIFF-verified 56M_5
triangle, whose sides measure 65.0, 72.6 and 72.6 m. It requires a simple closed
outline with three corners, sides in the 56–84 m range, and a shortest/longest
side ratio of at least 0.8. Extra edge vertices are allowed, with at most 1 m of
outline deviation and 3% area difference from the fitted triangle. Rotation,
winding direction, starting vertex, source ID and elevation do not determine
classification. Accepted triangles become polygon centroids and retain `VALUE`.
The exported spot properties include the measured side lengths and fit errors.

`--triangle-side 70` sets the nominal side length; its size and outline
tolerances scale with it. `--triangle-side 0` disables automatic triangle
conversion. Classification remains a heuristic: real contours matching a
symbol's shape and size can pass too.

Some source symbols are split into arcs at contour intersections. Before
classification, the script joins short fragments with exactly matching endpoints
and identical elevations. It accepts a reconstructed ring only if it passes the
same circle or triangle tests. It does not snap endpoints or fill missing edges.
All contributing source attributes and component hashes are saved in
`source_fragments`, with `classification_method` set to
`geometry_heuristic_joined_fragments`.

Optional `--spot-review` corrections match the normalized source geometry hash, source
ID, and elevation. They convert only the specified closed symbols to polygon
centroids. A changed or missing source geometry fails validation. Review files
are sheet-specific and do not change source elevation values.

The selected sheet grid uses **WGS84 (EPSG:4326)**, as agreed after the datum
comparison. DEM rasters use the corresponding WGS84 UTM zone (EPSG:32643 for
53A_16). Everest/Kalianpur boundary comparisons remain diagnostic evidence and
do not change the processing grid.

Contours are clipped in WGS84 before projection to avoid losing lines on sheet
edges. Original vertices are retained and long segments are sampled at half the
output pixel spacing. Missing, nonfinite, and conflicting coincident elevations
are rejected; zero is a valid elevation. Exact duplicate geometries at the same
elevation are removed.
Zero-length source lines are counted and excluded from interpolation. Their
original records remain in the cached source extract.

By default, `--edge-sliver-m 2` excludes an uncut line part only when more than
half its length is outside the selected sheet and its entire retained portion
lies within the 2 m inward boundary band. Lines entirely inside that band are
kept, as are lines reaching farther into the sheet. Point symbols are classified
first and are not subject to this filter. `--edge-sliver-m 0` restores strict
polygon clipping without sliver removal. This prevents neighbouring-sheet
fragments from supplying contradictory constraints to a single sheet; it does
**not** repair cross-sheet height discontinuities for a mosaic.

Both scripts now default to `--form-line-interval auto`. The detector joins
exact same-height fragments and finds open branches ending at an original vertex
of a continuing contour. It assumes one regular interval per selected sheet and
infers it from that sheet's contour heights, independently of the suspect
branch's height. Different sheets can have different intervals. There is no
universal 20 m assumption or requirement that the elevation series start at zero.

Reference evidence uses contour length inside the selected sheet, excluding
point symbols and suspect branches. A reference level needs 10% of the strongest
level's length. Adjacent height gaps vote in proportion to the shorter level's
length. The leading gap needs a majority of that vote; at least three levels and
90% of reference length must fit the inferred regular elevation series. Sparse
or inconsistent evidence leaves the interval unknown. These thresholds are
heuristic, not statistical confidence.

Only a branch half the inferred interval below every regular-series parent is excluded.
Uncertain branches, other off-interval contours, endpoint-to-endpoint height
changes, circles and triangles remain included. The reference-length thresholds
never remove a contour by themselves. Detection uses uncut vectors after symbol
classification and requires no toposheet raster. This recognises the calibrated
47J_9 junction pattern, not all possible form-lines or cartographic merges.
`--form-line-interval 20` supplies a known interval for the sheet;
`--form-line-interval none` disables detection.

Complete excluded line parts, source heights and geometry hashes are saved to
`excluded_contours.geojsonl` with the reason and supporting measurements. These
filters do not change heights. Reports explicitly record `seam_reconciled: false`.

Both scripts default to `--height-recovery auto` before contact handling and
interpolation. Two differently valued chains must join exactly at both ends to
form an isolated, simple closed loop. Enclosing and inner regular contours can
select one existing height using the inferred sheet interval. If those references
are missing or incomplete, a fallback can correct the anomalous chain only when
its complete connected length is at most **10 m and 10% of the loop perimeter**.
It identifies the regular-series height already joined at both ends, provided any
available enclosing/inner references agree. Crossings, other contacts, branching,
contradictory references and proximity to the sheet boundary prevent recovery.
No raster, reference DEM or manual case list is used.

`--height-repair-max-length-m` and `--height-repair-max-fraction` set the fallback
limits; zero for either disables the fallback. These caps do not restrict the
separate recovery based on complete nesting evidence. `--height-recovery none`
disables both rules. With the default Imagico engine, these rules select suspect
fragments whose heights are **unassigned**; their geometry is retained with an
unknown-contour ID. The preliminary surface estimates one constant height per
ID, which becomes a fixed constraint for the remaining stages. It is not snapped
to the sheet interval. The linear baseline still applies the proposed vector
repair directly. Original heights and selection evidence are retained in
`height_recovery.geojsonl` and `unknown_height_contours.geojsonl`; successful
Imagico runs write the estimated heights to `recovered_contours.geojsonl`.
The original classified `contours.geojsonl` and source extracts retain their
heights. All height points remain DEM inputs.

Regular contours can touch because their separation cannot be drawn at map
scale. The DEM builder recognises isolated tangencies and short lenses between
two continuing contours. Lens length along each contour must be at most twice
the relaxation radius. The default radius is half the output grid spacing
(15 m at 30 m resolution); `--contour-contact-radius-m` overrides it, and zero
disables handling. Only interpolation sections inside those small regions are
relaxed. Both source elevations and the complete classified contour export are
retained; all height points remain inputs, including points inside a contact
region. The operation does not recover the true ground separation.

Endpoint height changes, T-junctions, isolated crossings, shared segments and
contacts near the outer boundary remain strict conflicts. A connected contour
cannot be removed completely by this operation. The contact audit and relaxed
sections are exported separately, including when other conflicts prevent a DEM.

The default engine rasterizes contours and all height points onto a square UTM
grid. Where different contours occupy one pixel, the closest sampled line to
its centre wins; these collisions are counted. A point's height takes priority
in its pixel, and conflicting point heights in one pixel require a finer grid.
Known and estimated unknown constraints are restored exactly after refinement.
The directional search uses 32 rays extending up to 120 pixels. Cells with no
initial ray support remain nodata, and output is clipped to the sheet. The
linear baseline additionally limits output to the sample convex hull. Neither
method establishes that gaps between source constraints contain valid terrain.
An existing nonempty output directory is refused; use a new directory per run.

Outputs:

| File | Purpose |
| --- | --- |
| `dem.tif` | Float32 UTM elevation grid with -9999 nodata; includes spot candidates |
| `dem_contours_only.tif` | Validation run with point heights held out; never the final DEM |
| `imagico_validation_report.json` | Validation errors and unknown-height estimates before points are restored |
| `imagico_report.json` | Final numerical stages, raster collisions, unknown estimates, timing and implementation provenance |
| `interpolation_inputs.npz` | Final known-height raster and unknown-contour ID raster |
| `dem_water_flattened.tif` | Optional variant with eligible water polygons flattened |
| `constraint_distance.tif` | Distance in metres to the nearest sampled elevation constraint |
| `contours.geojsonl` | Classified, clipped contours in WGS84 with original heights |
| `recovered_contours.geojsonl` | Working contour copy used for interpolation, with audited height repairs |
| `unknown_height_contours.geojsonl` | Suspect fragments with null `VALUE`, original height, geometry and selection evidence |
| `height_recovery.geojsonl` | Candidate loops, decisions, original/recovered heights, reference evidence and fragment sizes |
| `spots.geojsonl` | Both circular and triangular symbols converted to WGS84 Points |
| `triangle_points.geojsonl` | The triangle subset, with the same properties as `spots.geojsonl` |
| `excluded_contours.geojsonl` | Complete source line parts excluded as outside slivers or form-line candidates, with heights and reasons |
| `form_line_candidates.geojsonl` | Branch candidates, including retained uncertain cases, with sheet interval and source provenance |
| `contour_contacts.geojsonl` | Cartographic contacts whose local geometry was relaxed, with both heights and the radius |
| `relaxed_contour_sections.geojsonl` | Original line sections relaxed for interpolation, retaining heights and provenance |
| `classification.json` | Input hashes, classification counts, height-recovery limits and decisions, including on failed DEM runs |
| `conflicting_constraints.geojsonl` | Coincident sample locations and their differing elevations, when interpolation fails |
| `report.json` | Input hashes, counts, grid metadata, coverage, residuals, and water decisions |

Classification vectors are written before interpolation, so the height-point
dataset survives a DEM failure. A directory containing these files without
`report.json` is not a completed DEM run. Conflicting elevations left after the
audited recovery rules are reported and prevent interpolation; they are not
averaged.

The point layers retain source attributes (`id`, `OBJECT_ID`, `TYPE`, `VALUE`,
etc.) and add:

| Property | Meaning |
| --- | --- |
| `symbol_type` | `circle` or `triangle`; native source Points use `point` |
| `elevation_m` | Numeric source elevation, without rounding or map-label corrections |
| `source_id`, `source_sheet`, `source_part_index` | Source feature, selected sheet and zero-based multipart component |
| `source_geometry_sha256` | Hash of the normalized original component in WGS84 |
| `source_fragments` | For reconstructed symbols, every original component's attributes, source ID, part index and geometry hash; the top-level hash describes the reconstructed outline |
| `classification_method` | `geometry_heuristic`, `geometry_heuristic_joined_fragments`, `source_geometry`, or `manual_review` |
| `symbol_radius_m`, `symbol_roundness` | Circle measurements |
| `symbol_side_min_m`, `symbol_side_max_m` | Triangle side measurements |
| `symbol_edge_error_m`, `symbol_area_error_fraction` | Triangle fit measurements |

`spot_classification` retains the more detailed detection label used by earlier
runs. Manually supplied symbols whose shape fails both automatic tests use
`symbol_type: "unknown"`. Symbol measurements are in the sheet's UTM projection.
The [SOI legend](https://storage.googleapis.com/soi_data/legend.webp) identifies
triangles as triangulated height stations. Some circular candidates correspond
to benchmark labels in the calibration maps; the automatic export deliberately
records circle versus triangle without inferring benchmark subtypes or peaks.

Water polygons are unioned before processing. Enclosed polygons with fully
supported shorelines and at most 5 m of modelled shoreline variation are flattened
to the median modelled shoreline elevation. Polygons crossing the sheet boundary,
unsupported polygons, and variable shorelines are skipped and recorded. This is
an explicit approximation, not a measured water level or river-conditioning step.
The unflattened DEM is always retained for comparison.
These CartoDEM water polygons come from a different source. They provide a
shared smoothing mask for consistency with the other DEM dataset; they are not
evidence for accepting or rejecting height points or changing their elevations.

## Toposheet calibration: 56M_5

The [GeoTIFF listing](https://github.com/ramSeraph/india_topo_maps/releases/download/50k-osm-georef/listing_files.csv)
locates this sheet in the `50k-osm-georef-extra4` release. Use its listed URL:

```sh
curl -fsSL \
  https://github.com/ramSeraph/india_topo_maps/releases/download/50k-osm-georef-extra4/56M_5.tif \
  -o data/contours/56M_5/toposheet.tif

uv run audit_spot_heights.py data/contours/56M_5/contours.geojsonl \
  data/contours/56M_5/toposheet.tif --sheet 56M_5 \
  --output-dir data/contours/56M_5/calibration-auto
```

The audit saves raw crops, crosshair crops, contact sheets, shape measurements,
and input hashes. `--context-m 1100` with a new output directory produces wider
crops for labels extending beyond the default 600 m window. Blue crosshairs mark
vector centres; they do not imply a measured registration accuracy.

Manual review of all 84 circle candidates and 24 small rejected geometries is
recorded in [calibration/56M_5.review.csv](56M_5.review.csv):

- 63 circles have approximately 8.75 m radius and 21 approximately 11.25 m.
  All have roundness between 0.999036 and 0.999558 and correspond to height
  annotations in the map. Larger circles include benchmarks marked `BM`.
- 81 circle heights were read and agree with `VALUE` within 0.001 m. Three
  benchmark labels were only partly visible in the reviewed crops; their
  numeric heights remain unverified. The review distinguishes these cases.
- Feature `CONTOUR_Crv.4519423` is a triangular symbol labelled 177 m. The
  automatic triangle detector now finds it without a reviewed correction.
  [calibration/56M_5.spots.json](56M_5.spots.json) retains the earlier
  manual correction and GeoTIFF provenance for comparison. Its crop number
  refers to the original `calibration-final` audit, before automatic triangles
  were grouped with the detected symbols.
- The other 23 reviewed controls are contour lines or fragments. This is a
  targeted review, not an exhaustive search of the raster for missing heights.

Automatic classification agrees with all **108 reviewed features**: 84 circles,
one triangle and 23 contours. The resulting DEM has **739 contour parts and
85 point heights**, with no review file supplied. The original 84-point baseline
is in `data/contour_dem/56M_5`, the manually corrected comparison is in
`data/contour_dem/56M_5-calibrated`, and the automatic run is in
`data/contour_dem/56M_5-auto-symbols`. The GeoTIFF is a rendered map from the same
survey family: it validates symbol interpretation, not independent terrain
accuracy.

## Cross-sheet calibration: 56M_9

The unchanged detector finds 92 circles and 12 triangles on 56M_9, leaving
598 contour parts. All 12 triangles and their printed elevations were verified
against the [56M_9 GeoTIFF](https://github.com/ramSeraph/india_topo_maps/releases/download/50k-osm-georef/56M_9.tif).
Of the 92 circles, 88 printed heights agree, two labels disagree with the source
values, and two map associations remain uncertain. Twenty-four small closed
contours were also reviewed and correctly retained as lines. Full observations
are in [the review CSV](56M_9.review.csv) and
[validation record](56M_9.validation.json).

The 178 m circle within the mapped lake has the same approximately 8.75 m radius
as the other symbols and remains included. All source elevations are retained;
the calibration notes do not become runtime map dependencies or corrections.
This candidate review checks interpretation and label agreement, not exhaustive
detection recall or independent elevation accuracy.

The current exports with symbol properties are in
`data/contour_dem/56M_5-symbols` (84 circles, one triangle) and
`data/contour_dem/56M_9-symbols` (92 circles, 12 triangles). Each contains both
`spots.geojsonl` and `triangle_points.geojsonl`.

## Five additional sheets across India

The same symbol size and shape thresholds were applied to five more sheets.
The [results record](cross_country_5.results.json) links each run,
review sample and coverage exclusion. All runs use the released GeoParquet
extracts; the GeoTIFFs are calibration inputs only. No water mask was applied.

| Sheet | Region | Circles | Triangles | DEM outcome |
| --- | --- | ---: | ---: | --- |
| 53A_16 | Himachal Pradesh | 155 | 0 | 37 conflicting sample locations |
| 45F_4 | Rajasthan | 140 | 10 | Created; 98.17% sheet coverage |
| 73D_15 | Odisha | 101 | 21 | Created; 99.65% sheet coverage |
| 83E_11 | Arunachal Pradesh | 84 | 14 | 10 conflicting sample locations |
| 47J_9 | Maharashtra | 95 | 15 | 46 conflicting sample locations |

All **635 height points**, including **60 triangles**, are saved even where
interpolation fails. The completed DEMs include every detected height point.
The contour-only holdout RMSE is 5.58 m for Rajasthan and 5.18 m for Odisha;
these measure internal agreement, not independent terrain accuracy.

Himachal exposed circles split into fragments at contour crossings. Joining
exact same-elevation endpoints recovered 19 additional circles, all confirmed
as height points in the raster, and reduced conflicting sample locations from
77 to 37. Eighteen of those 19 printed heights agree; one source value is
1419 m where the map reads 1423 m. The source value is retained.

The Himachal, Odisha and Maharashtra reviews mostly support the symbol
interpretation, with a few differing or unreadable labels and one uncertain
map association. Rajasthan has many missing or differing raster annotations;
Arunachal's raster omits the reviewed contours and height labels. Those two
rasters cannot establish detector accuracy. These are sampled reviews, not
an exhaustive search for missed symbols.

The remaining Himachal conflicts are present in the original released
geometries, before our sheet clipping. A candidate Kalianpur 1962 to WGS84
transformation places 18 of 37 within 25 m of the sheet or internal 7.5-minute
grid, compared with 10 for the unshifted grid. This supports investigating
datum-related seams but does not establish the original cutting datum.
Interior conflicts remain, including very short lines with incompatible
elevations. No coordinate transformation or automatic elevation correction has
been applied. Detailed measurements are in
`data/contours/53A_16/datum-boundary-analysis.json`.

### Contour-line investigation

Tracing the remaining conflicts to the uncut released geometries found three
different patterns. All 93 locations match actual source intersections; none
are introduced solely by our clipping or millimetre sample rounding.

| Sheet | Observed source geometry |
| --- | --- |
| 53A_16 | 29 two-line joins with changing heights; 8 four-way endpoint junctions |
| 83E_11 | 10 two-line joins, all within 2 m of the nominal WGS84 outer sheet edge |
| 47J_9 | 46 T-junctions involving 36 lower-valued branches, each ending against a contour 10 m higher |

Himachal includes a 2.74 m segment labelled 400 m between two long 1,000 m
sections, a 6.63 m segment labelled 1,240 m closing an otherwise 880 m loop,
and a 5.95 m segment labelled 1,360 m closing an otherwise 1,060 m loop.
These are strong candidates for incorrect fragment attributes, not automatic
corrections. Eight conflicts occur at four small lenses where differently
labelled lines meet twice. The cutting-grid hypothesis still fits some other
Himachal joins, including internal quarter-sheet seams. Reinterpreting the
outer index as Kalianpur leaves 29 of the 37 existing locations inside; this
is not a complete re-extraction or a count of conflicts in a repaired dataset.
Arunachal's joins track the unshifted outer grid more closely than either
tested Kalianpur grid. Maharashtra's conflicts lie at least 215 m from all
three tested outer/quarter grids and all remain inside either shifted sheet.

The 47J_9 raster reveals that the lower-valued branches are drawn as dashed
relief lines. For example, C004 has a solid contour visibly labelled 700 m,
and a dashed branch meeting it whose released `VALUE` is 690 m. C007 shows
a dashed spur with source value 730 m joining a solid 740 m contour. The
vectors align with these map features; a sheet-grid shift does not explain
these junctions. All 46 locations have saved map crops; some exact joins are
partly covered by text, roads, buildings or other symbols.

The [SOI Handbook of Topography, Chapter V, para 36 (PDF page 55)](https://surveyofindia.gov.in/documents/soichapter-v.pdf#page=55)
describes fine broken form-lines as relief detail that can occupy an
intermediate altitude rather than a prescribed contour level. The
[Chapter VI feature catalogue (PDF page 93)](https://www.surveyofindia.gov.in/documents/soichapter-vi.pdf#page=93)
lists broken form-lines/sub-features separately from regular contours and
cliffs. The raster pattern supports a form-line interpretation, but the
released attributes provide `TYPE=1` for both classes and no such feature
code. A rocky outcrop may be the landform depicted; that does not establish
that the dashed line's numeric value is an exact constant-height constraint.
These findings do not rule out cliff-related contacts elsewhere.

The initial investigation did not alter production constraints. The 36
implicated Maharashtra branches are preserved, with their original properties
and geometry, in `calibration/47J_9.form-line-candidates.geojsonl`. The subsequent
optional detector described above identifies exactly those 36 source parts,
with no additional or missing matches on 47J_9. The combined evidence record is
`calibration/contour-conflicts.findings.json`.

### Single-sheet filtering and cross-sheet seams

The user chose to exclude identified form-lines from DEM constraints and permit
clear neighbouring-sheet slivers to be ignored for individual sheets. New runs
with these policies retain all circular and triangular height points:

| Sheet | Applied filters | Result |
| --- | --- | --- |
| 83E_11 | 2 m edge band; 739 mostly-outside source parts excluded | DEM created, 99.9979% coverage, all 98 height points included |
| 47J_9 | 2 m edge band plus 20 m form-line interval; 44 outside parts and 36 form-line candidates excluded | DEM created, 99.9548% coverage, all 110 height points included |
| 53A_16 | 2 m edge band | 34 conflicting sample locations remain, down from 37; DEM still refused |

The new output directories are `data/contour_dem/83E_11-edge-filtered`,
`data/contour_dem/47J_9-form-lines-excluded`, and
`data/contour_dem/53A_16-edge-filtered`. These are single-sheet results, not a
validated seamless mosaic. The older runs above remain unchanged as baselines.

The initial local automatic-interval check covered seven cached sheets; see
`calibration/form-line-intervals.findings.json` and the
`data/contour_dem/*-auto-interval-audit` directories. On 47J_9 it excludes 34 of
the original 36 reviewed candidates plus one same-height continuation fragment
(`2319362`), also confirmed dashed in a separate raster review. Two reviewed
branches (`2319356`, `2319418`) have inconsistent local spacing evidence and
remain constraints, leaving two conflicting locations. No form-line exclusions
are added on the other six sheets. All 824 point coordinates, heights, symbol
types and source geometry hashes match their previous exports. These are
classification and sample-conflict checks retained as historical results.

The current sheet-level interval check is recorded in
`calibration/sheet-intervals-and-contacts.findings.json`, with outputs under
`data/contour_dem/*-sheet-interval-audit`. It infers 10 m for 45F_4 and 20 m for
the other six sheets. On 47J_9 all 36 reviewed candidates and the connected
continuation fragment are excluded, resolving the two previously uncertain
cases. No form-line exclusions are added on the other six sheets. Cartographic
contact handling relaxes eight locations forming four small lenses on 53A_16,
leaving 26 other sample conflicts. The other six sheets have no sample conflicts.
All 824 height points still match the previous exports. The 53A_16 review queue
is `calibration/53A_16.remaining-conflicts.review.json`.

The initial nesting-only height-recovery check is retained in
`data/contour_dem/53A_16-height-recovery-audit-checked` and
`calibration/height-recovery.findings.json`. It recovered C001/C002 to 1,260 m,
C003/C004 to 1,360 m and C019/C020 to 1,880 m, changing four source parts and
reducing 26 conflicts to 20. Multiple inner loops were allowed when their
heights agreed.

The current production rerun includes the size-limited fallback described above:

```sh
uv run contours_to_dem.py data/contours/53A_16/contours.geojsonl \
  --sheet 53A_16 --output-dir data/contour_dem/53A_16-height-autocorrect
```

It additionally recovers three fragments, removing six conflict locations:

| Cases | Original → recovered height | Connected fragment length | Fraction of loop |
| --- | --- | --- | --- |
| C005/C006 | 1,240 → 880 m | 6.63 m | 5.13% |
| C013/C014 | 400 → 1,000 m | 2.74 m | 0.022% |
| C025/C026 | 1,360 → 1,060 m | 5.95 m | 0.743% |

Together the two rules recover six loops across seven source parts. **14
conflict locations remain**, so this run still refuses to create a DEM. All
155 height points and eight cartographic contact decisions match the baseline.
The other six cached sheets produce no recovery proposals or height changes.
The 109 passing checks cover recovery, size limits, fragmented chains,
interpolation, extraction, form-line inference and cartographic contacts;
findings are in `calibration/height-autocorrect.findings.json`.

The local GLO-30 check in `data/comparisons/53A_16-C005-glo30/report.json` found
868–870 m at the C005/C006 endpoints and a maximum of 888 m within 100 m of
C005, supporting the 880 m interpretation. Its roughly 30 m sampling cannot
exclude an extremely narrow spike by itself. This was calibration evidence and
is not an input to the recovery rule.

`audit_height_recovery.py` remains available for standalone diagnostics using
the same recovery function. Subsequent visual review uses only
`data/contour_dem/53A_16-height-autocorrect/interior_review.json`, containing six
interior conflict locations. Eight outer-boundary cases or candidates are kept
separately in `boundary_review.json`; all 14 remain in `unresolved_review.json`
and in the numerical conflict audit. This review separation does not resolve or
remove any DEM constraints. Recovered locations stay in the review history but
are excluded from the active queue.
`plot_conflict_neighbourhood.py` provides height-coloured raster overlays and
contour-only views.

The C011/C012 trace in
`data/contours/53A_16/line-conflict-audit/C011-recovery-trace.json` exposes a
limitation of both current recovery rules: the 19.99 m, 740 m fragment joins
two different open 760 m chains, so it is not a closed-loop candidate and never
reaches the nesting checks. Its length also exceeds the 10 m fallback limit.
The automatic selector still leaves it unresolved. An explicit vector-only
calibration with `audit_imagico_case.py` unassigned this fragment's height and
estimated **758.84 m at 5 m grid spacing** and **759.03 m at 10 m**, close to the
adjoining 760 m contours. Both final runs include the one local height point.
At 30 m, the fragment shares its only pixel with a known 760 m contour: the
reported 760 m estimate adds no independent constraint pixels, so it is not
evidence of greater accuracy. These experiments do not change the automatic
selector or its 10 m fallback cap, and do not edit the source data.

```sh
uv run audit_imagico_case.py data/contour_dem/53A_16-height-autocorrect \
  --sheet 53A_16 --unknown-source CONTOUR_Crv.2490643 \
  --resolution 5 --output-dir data/contour_dem/53A_16-C011-imagico-5m
```

Use a new output directory to repeat this command. This diagnostic selects
source IDs explicitly; the production pipeline has no hardcoded case list and
requires no toposheet raster. Full-sheet `53A_16-imagico-preflight` still stops
at 14 unresolved known-height conflicts after unassigning seven source parts
in the six supported recovery cycles. It does not produce a DEM for that sheet.

The 53A_16 crop collection covers the original 37 conflicts. Its review CSV now
marks C018, C029 and C030 as removed by the edge filter, eight contact locations
as locally relaxed, 12 locations as height-recovered and 14 as remaining.
The 1,080/1,100 m C023 junction is 5,591 m from the outer boundary and
4,642 m from the nearest internal quarter boundary, so it is an interior case.
Detailed crops show both solid-contour attribute changes and solid contours
meeting; the user identified C023 as regular contours touching at map scale.

`audit_contour_seams.py` compares the heights on both sides of exact original
endpoint joins within 2 m of the nominal outer sheet boundary, including matching
heights as controls. It uses the majority of each complete source line's length
to infer its side of the sheet; source sheet IDs are absent from the release.
The 83E_11 cached extract contains 829 such unambiguous paired locations:

| Edge / nominal neighbour | Compared joins | Agreement | Mismatches, selected sheet minus neighbour |
| --- | ---: | ---: | --- |
| West / 83E_7 | 332 | 330 | Two at +20 m |
| South / 83E_12 | 308 | 300 | Two at +20 m, two at -100 m, four at -10 m |
| North / 83E_10 | 189 | 189 | None among these paired endpoints |
| East / 83E_15 | 0 | — | Not assessed by this cached extract |

All ten mismatches involve a neighbouring line whose retained portion is confined
to the 2 m edge band. Removing those slivers clears the single-sheet conflicts,
but the original heights still disagree across the seam. The many agreeing joins
and differing mismatch sizes rule out a single constant height offset as a fix
for this sample. They do not establish which member of each mismatched pair has
the correct height, or whether a line was joined to the wrong contour.

53A_16 has 617 paired outer-edge locations in the cached extract; 613 agree and
four differ (-20 m once, +20 m twice, +40 m once). This outer-edge audit does not
cover its internal quarter-sheet seams or other interior conflicts.

Each audit saves `report.json`, `joins.geojsonl`, `source_lines.geojsonl`, and
`height-jumps.png` under `data/contours/SHEET/seam-audit`. These records preserve
the unresolved source disagreements for mosaic work. They cover only exact
paired endpoints available in the cached extract, not a complete survey of
adjacent sheets or nearby unmatched endpoints.

```sh
uv run audit_contour_seams.py data/contours/83E_11/contours.geojsonl \
  --sheet 83E_11 --output-dir data/contours/83E_11/seam-audit
uv run extract_sheet_dem.py --sheet 83E_11 --reuse-extracts --no-water \
  --output-dir data/contour_dem/83E_11-edge-filtered
uv run extract_sheet_dem.py --sheet 47J_9 --reuse-extracts --no-water \
  --form-line-interval 20 --output-dir data/contour_dem/47J_9-form-lines-excluded
```

Reproduce the line audit and raster crops using cached extracts:

```sh
uv run audit_contour_conflicts.py data/contours/47J_9/contours.geojsonl \
  data/contour_dem/47J_9-cross-country/contours.geojsonl \
  data/contour_dem/47J_9-cross-country/conflicting_constraints.geojsonl \
  --sheet 47J_9 --output-dir data/contours/47J_9/line-conflict-audit
uv run audit_conflicts_toposheet.py \
  data/contours/47J_9/line-conflict-audit/report.json \
  data/contours/47J_9/toposheet.tif \
  --output-dir data/contours/47J_9/line-conflict-audit/toposheet-review
```

`line-conflict-audit` contains source IDs, exact-contact types, lengths,
endpoint neighbours, candidate-grid distances, source-vector plots and an
overview. `toposheet-review` contains native raster crops, paired raw/overlay
images, six contact sheets and a visual-review CSV. These audits cover the
reported sample conflicts, not every possible unsampled contour crossing.

The initial coverage probes 47F_9, 83B_10, 73H_3 and 56N_14 yielded empty
extracts or only boundary fragments; they were replaced and remain recorded
as coverage gaps. Partition bounding boxes alone do not establish sheet coverage.

## Validation and interpretation

Compare a completed sheet with standard Copernicus GLO-30:

```sh
uv run compare_dem.py data/contour_dem/73D_15-cross-country/dem.tif \
  --spots data/contour_dem/73D_15-cross-country/spots.geojsonl \
  --constraint-distance data/contour_dem/73D_15-cross-country/constraint_distance.tif \
  --output-dir data/comparisons/73D_15-glo30-native-grid
```

This reads only the required windows of the AWS Copernicus DEM 2021 COGs,
preserves their native pixel origin, and bilinearly resamples onto the exact SOI
grid. It saves the native reference window, aligned reference, difference raster,
point comparison CSV, plot and a report with input hashes. Cached reference
windows are checked before reuse. Source height points are sampled at their
exact coordinates in the native reference window, separately from the DEM grid.
No vertical conversion, fitted bias correction or horizontal shift is applied.

On 73D_15, the final comparison covers 800,173 common 30 m pixels (720.16 km²).
All differences below are **SOI minus GLO-30**, in metres:

| Samples | Mean | Median | MAE | RMSE |
| --- | ---: | ---: | ---: | ---: |
| All common pixels | -0.12 | +0.14 | 5.62 | 7.23 |
| All 122 original height points | +4.26 | +2.91 | 5.62 | 8.17 |
| 31 points on GLO-30 slopes below 5° | +0.31 | +0.28 | 3.47 | 5.43 |

The absence of a large overall offset is consistent with heights referenced
near mean sea level, but it does not establish that SOI uses EGM2008. The
[SOI geodesists' description of topographic contours](https://www.fig.net/resources/proceedings/fig_proceedings/fig2018/papers/ts06e/TS06E_singh_srivastava_9496_abs.pdf)
states that their heights are MSL based. The precise Indian vertical datum
realization of these released `VALUE` attributes remains unverified.
[Copernicus specifies](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM)
WGS84 horizontally and EGM2008 vertically. WGS84 coordinates alone do not say
that elevation values are ellipsoid heights.

Copernicus is a surface model that includes vegetation and buildings. Local
differences also reflect interpolation, registration, source epochs and surface
detail; these statistics are not an independent estimate of SOI terrain accuracy.
The initial unaligned-crop experiment is retained in
`data/comparisons/73D_15-glo30` and marked superseded; use the `native-grid` result.

## Python Imagico validation

Two cached sheets were rerun with the default Python engine at 30 m, without
toposheet inputs or water smoothing:

| Sheet | Final circle / triangle points | Contour-only holdout | Final supported sheet cells |
| --- | --- | --- | --- |
| 56M_9 | 92 / 12 | 104 predictions; 4.55 m RMSE | 100% |
| 45F_4 | 140 / 10 | 140 predictions; 5.73 m RMSE; 10 unavailable | 97.57% |

All 254 points are included in the final runs, including the ten without a
validation prediction. Sampling their exported raster cells reproduces their
source heights within Float32 precision (maximum error 0.0000062 m). Final
numerical interpolation took approximately 28 and 25 seconds respectively,
plus a separate validation pass and vector preparation; these are individual
local runs, not controlled benchmarks. Results are in
`data/contour_dem/56M_9-imagico`, `data/contour_dem/45F_4-imagico` and
`calibration/imagico.findings.json`. This is initial calibration, not proof of
accuracy across all terrain.

```sh
uv run test_contours_to_dem.py
uv run test_contour_form_lines.py
uv run test_contour_contacts.py
uv run test_contour_height_recovery.py
uv run test_extract_sheet_dem.py
uv run test_compare_dem.py
uv run test_imagico_interpolate.py
uv run test_imagico_dem.py
```

Tests cover a known planar surface through projection and raster export, valid
zero elevations, missing/conflicting elevations, open and chord-closed arcs
versus circular symbols, triangle shape and scale rejection, rotated and
densified triangles, the real 177 m triangle without a source-ID exception,
classification before multipart clipping, reviewed corrections, point-property
preservation, separate triangle export, and water exclusions. Pipeline tests
exercise extraction through an iomaps stub, cache provenance checks, empty water
extracts, and failed downloads; the actual sheet runs use cached iomaps extracts.
Additional tests check that boundary filtering preserves deeper crossings,
inside lines and height symbols, and that form-line filtering preserves ordinary
half-interval lines, endpoint-to-endpoint joins and height points.

The report evaluates a **contour-only** interpolation at the held-out spot
locations, then reruns the **complete pipeline with every point included**.
Unknown contour heights are recomputed using those points; the validation
estimates are not reused. Both circle and triangle heights are ordinary final
constraints. This measures internal source agreement, not independent survey
accuracy. Contour midpoint residuals measure
reproduction of the input lines between samples; they are not a ground-truth
accuracy estimate either. A regression test with points above a flat contour
surface verifies that both circle and triangle elevations affect the final DEM
after the contour-only holdout check.

30 m is grid spacing, not a claimed horizontal or vertical accuracy. The Python
Imagico adaptation still needs broader terrain calibration; it is not
hydrologically conditioned. Rasterization can merge closely spaced constraints,
and exact preservation applies at constraint pixels, not every subpixel vector
coordinate. Source vertical datum, spot omissions, data gaps and joins between
sheets remain separate validation concerns. Historical TIN measurements above
do not describe the new engine unless a run is explicitly labelled Imagico.
