# WRIS river selection and constrained joins

Historical WRIS record. Production subsequently returned to GLO-30 river class 3
to cover the full CartoDEM footprint outside India as well. See
[GLO30_INTEGRATION.md](GLO30_INTEGRATION.md). The shared-junction fitting method
developed here is retained; the WRIS measurements below remain historical.

Earlier method revision: 24 October 2026. WRIS replaced GLO-30 class 3 as the river
footprint. GLO-30 continues to supply only the ocean boundary. River heights
come from EGM2008-corrected CartoDEM. The two-sheet pilot is a method-validation
run; it does not establish nationwide coverage or justify publication by itself.

**Correction:** `pilot-wris-v2` used an overly strict per-pixel source-height
constraint that was not the agreed shared-boundary model. Its 6.2% coverage is
a rejected experiment, not completion of the river method. The current `joint`
method and its separate `pilot-wris-v3` validation are described below. Older
source snapshots preserve their original notes for reproducibility.

## Source and feature classes

Use the full-resolution GeoParquet from the
[river release](https://github.com/ramSeraph/indian_water_features/releases/tag/rivers),
not simplified vector tiles. The iomaps catalogue route is
`/river-polygons/wris/`. `wris.lock.json` pins the download URL, byte size and
SHA-256. `build.py download` fetches/verifies it; the river geometry stage also
fetches it if needed. The verified source contains 40,898 nonempty polygons.

The [national ArcGIS layer](https://arc.indiawris.gov.in/server/rest/services/SubInfoSysLCC/WaterBodies/MapServer/9?f=pjson)
uses a single symbol and has no coded-value domain. The unique-value renderer of
[WRIS's Ganga Basin river-polygon layer](https://arc.indiawris.gov.in/server/rest/services/NWIC/River_GangaBasin/MapServer/11?f=pjson)
explicitly labels the same `ripcode` field:

| Value | Label | Features in pinned Parquet | Treatment |
| --- | --- | ---: | --- |
| `01` | River | 9,837 | Candidate river footprint |
| single space | River | 91 | Candidate river footprint |
| null | River | 112 | Candidate river footprint |
| `02` | Sand | 24,537 | Excluded, source heights retained |
| `03` | Island | 6,321 | Excluded, source heights retained |

Null and blank are included because the renderer explicitly labels them River;
this is not a general rule that missing attributes mean water. Unknown codes
fail. The live service has two additional blank-code rows, OBJECTIDs 281/282,
with empty geometries, explaining its 40,900-row count. All eight source
attribute fields survive in our Parquet. `bacode` is entirely null and `k50` is
almost entirely null; neither supplies direction information.

The saved field/legend evidence and three-area visual comparison are in
[`diagnostics/river-mask-comparison`](../diagnostics/river-mask-comparison/README.md).
Coastal mask agreement is much greater than at Patna or Guwahati. These are
agreement checks, not independent accuracy assessments or water-date checks.

## Raster masks and bounded reads

`wris.py` builds a disk-backed SQLite R-tree from batched WKB reads, preserving
the original polygon geometry and interior rings. It checks geometry type,
bounds, classes and unique source IDs. The cached index and per-sheet rasters
have input/code/output hashes. Candidate geometries are queried by sheet extent.
The pinned source has 45 invalid polygons. Self-touching linework may be
normalized with GEOS `make_valid(method='linework')`, retaining polygonal parts,
only when the absolute area change is less than one native pixel and the
original and normalized geometries rasterize identically on that sheet.
Any changed pixel membership fails. Per-sheet reports record source feature IDs,
invalidity reasons, original/repaired WKB hashes and the area comparison. The
original download and index remain unchanged. This does not authorize general
boundary repair or simplification.

Rasterization uses the exact native CartoDEM grid and pixel-centre membership
(`all_touched=False`). There is no buffer, resampling or simplification. River
classes are unioned; Sand and Island override any overlapping River geometry.
GLO-30 ocean then overrides river membership. GLO-30 river and lake classes do
not add to or subtract from WRIS River. Lake overlap counts are diagnostics only.
Outside WRIS River, existing CartoDEM elevations are retained. No separate lake
flattening is performed.

The internal class raster is 0 for other non-ocean pixels, 1 for GLO ocean, 3 for
WRIS River and 255 for unavailable classification. These internal numbers are
not WRIS `ripcode` values. In particular, WRIS `03` means Island, while internal
raster class 3 means the selected river footprint.

Per-sheet geometry uses a 512-pixel halo. Shared global pixel IDs assemble the
river network across sheet boundaries. The existing GRASS/RivGraph evidence,
support checks, shared junction stations and downstream constraints still apply.
The [national runner](NATIONAL_BUILD.md) uses overlapping contexts of at most
40 million cells, then fits complete connected components across all sheets.
The context boundaries do not partition the fitted river surfaces.

## Shared boundary elevations in the joint fit

Fit connected river sections together. A junction has one unknown elevation
used by every adjoining section; it is estimated from corrected CartoDEM bank
observations. It is not required to equal each noisy original water pixel.
This distinction matters: solver agreement at a shared model boundary can be
checked to numerical precision without claiming the source DEM is that accurate.

`joint.py` first fits the connected network with the independently supported
directions as constraints and curvature regularization on all reaches. This
allows neighbouring bank evidence and shared junctions to inform short or flat
sections whose direction was previously unresolved. For those sections, the
initial joint endpoint levels propose a direction. The second fit imposes
downstream monotonicity on every proposed directed section. Equal-level sections
share a level unknown; directed cycles forced by those equalities also collapse
to a common level. Every expanded Float32 reach is checked again.

Directions inferred this way are explicitly recorded as **joint profile
hypotheses**, separate from independent bank/GRASS evidence. This is an
experimental extension of the direction inference, not external confirmation of
flow. Fitting does not establish accuracy against surveyed water elevations.
Existing bank support, gap and adjustment checks remain in force.

The application covers complete modelled connected groups. If an unmodelled
native contact would still join a fitted part to an unresolved part, the group
remains unchanged; no independent fitted patch is pasted against raw heights.
All nodes and cross-sections use the same shared solution. This removes the
partial-treatment interface instead of imposing exact raw-pixel levels or an
arbitrary transition ramp. Final COG checks still inspect actual pixel gradients,
nonriver/unsupported preservation, nodata and sheet boundaries.

The default is `--join-method joint`. The rejected `source-pins` experiment is
available only as an explicit method option or from its saved source snapshot.

## Rejected experiment: exact source-pixel constraints

The historical v2 GLO-mask pilot reverted groups when new steps exceeded 10 m.
That threshold was a coarse diagnostic, not an acceptable artificial seam size.
The WRIS v2 experiment imposed source-pixel constraints in the profile fit.

For every pair of cardinally adjacent valid river pixels, including sheet
borders, identify interfaces between a fitted section and retained source data,
and between independently fitted groups. At each such interface, pin the
**treated pixel's own pre-fit elevation**, using its interpolation from shared
profile stations:

```text
(1 - fraction) * h[left_station] + fraction * h[right_station]
    = original corrected CartoDEM height at that treated pixel
```

The solver permits 0.5 mm around that equality and validates at most 1 mm after
Float32 storage and raster interpolation. This is numerical tolerance, not a
claim about CartoDEM accuracy or a physical maximum river slope. Holding the
correction to zero at the treated boundary preserves the original difference
between the two adjacent pixel centres. It does not copy the neighbouring
pixel's value or make the original local slope flat.

The same constrained fit retains downstream monotonicity, shared junction
elevations, bank upper bounds and the existing support/adjustment gates. If
different boundary samples, bank caps or direction constraints cannot support a
coherent flat cross-section and downhill profile, keep the entire affected
profile group unchanged. Rescan interfaces after rejection and refit remaining
groups. No isolated pixel is reverted within a fitted cross-section.

For that experiment, `profiles.py` accepts linear `boundary_samples`; `boundaries.py` finds interfaces
and runs the constrained refits; `directions.py` records each rejected group.
`audit.py` independently reads final native COG pixels to check interfaces within
and across sheets, as well as source validity, unchanged nonriver/unsupported
pixels and stored downstream station order.

That experiment prevents added correction jumps at partial-treatment interfaces. It
does **not** repair existing source steps, prove mask accuracy, or establish that
every accepted elevation is correct. More conservative support can reduce the
treated area. Tidal reaches and inconsistent directions can remain unresolved.

## Reproduction

Install `../requirements.txt` and the recorded GRASS CLI. After preparing the
coastal EGM2008 baseline, run from `egm2008/`:

```sh
python -m unittest -v test_pipeline test_production test_cog \
  river_flattening.test_profiles river_flattening.test_geometry \
  river_flattening.test_production river_flattening.test_directions \
  river_flattening.test_boundaries river_flattening.test_wris river_flattening.test_joint

python -m river_flattening.production run \
  --workdir work --label pilot-wris-v3 --river-source wris --join-method joint \
  --sheets cdnd43v.tif cdnd43w.tif

python -m river_flattening.audit --workdir work --label pilot-wris-v3 \
  --crop 75.08 12.10 75.24 12.28 \
  --crop-output work/rivers/pilot-wris-v3/coastal-crop.npz
```

The default source is WRIS. `--river-source glo30` is available for historical
mask comparisons; it does not restore the old fitter. Reproduce older pilots
from their saved source snapshots. Use a new output label after recipe changes.
Final Float32 COGs retain `NoData=-255`; neither river fitting nor routing creates
valid output in original holes. TerrainRGB follows the existing quantization
and alpha-254 nodata conventions after the COG validation.

## Pilot status

The source-pin experiment treated 30,029 of 481,308 valid WRIS river pixels
(6.24%). It rejected 105 groups due to incompatible exact source-height targets.
Its native-pixel audit passed, but those checks did not make the method appropriate
or establish completion of the shared-boundary plan. Files remain under
`pilots/wris-v2/` as a rejected experiment.

The separate joint-fit pilot completed the two-sheet COG export and native audit:
106,545 of 481,308 valid WRIS river pixels were treated (22.14%). All nonriver,
unsupported and nodata pixels remained unchanged. There are no partial-treatment
interfaces within the written river footprint, including sheet borders. Of 2,072
applied reaches, 538 have independent direction evidence and 1,534 use a joint
profile hypothesis. This distinction is material; a consistent fitted profile
does not independently establish flow direction.

Remaining concerns include newly increased **within-profile** adjacent-pixel
steps up to 9.59 m, despite the removal of partial-treatment joins. A short
station interval and discrete cross-channel owner assignments can still produce
steep pixel changes. The largest source-height correction is 64.22 m at the
previously discussed coastal-mask artefact. These results require further gradient
and surface-projection review before nationwide publication. Passing the native
preservation checks is not a claim of a fully validated river surface.

Reports, exact fit source, profiles and the small coastal diagnostic are in
`pilots/wris-v3/`. See [PILOT_RESULTS.md](PILOT_RESULTS.md) for the detailed checks.
