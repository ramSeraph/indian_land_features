# River direction integration

Implementation started 2026-10-23. The first integrated run uses both complete
D43V and D43W sheets, the same domain as the earlier 17.0% pilot. This document
describes implemented code; benchmark results belong in `PILOT_RESULTS.md`.
Nationwide validation remains pending.

**24 October update:** WRIS now supplies river geometry, and constrained refits
replace the 10 m post-fit gate for new builds. See
[WRIS_INTEGRATION.md](WRIS_INTEGRATION.md). This document preserves the earlier
GLO-30 pilot rationale and results; its 10 m acceptance section is historical.

## Processing unit

Ocean mouths anchor the downstream ends of connected rivers. Directions can
propagate inland through network continuity, supported by bank heights and
terrain routing. Sheets are storage/output partitions, not successive upstream
steps. A disconnected river-mask fragment cannot inherit an outlet through a
gap merely because its sheet adjoins a coastal sheet. Lake class 2 stays outside
the edited river network and outside bank sampling.

The profile solver still operates on shared stations across sheet boundaries.
It fits each connected group of accepted reaches once. An ambiguous branch can
remain unchanged while other supported parts of its network are treated. A
native skeleton support array prevents interpolation over an unresolved branch
whose endpoint junctions happen to have fitted elevations.

## External code actually used

- GRASS CLI `r.watershed -sm`, disk mode with a 256 MB segment-cache limit,
  generates D8 drainage and accumulation from corrected CartoDEM. The first
  integration uses GRASS **7.8.2**, Ubuntu package **7.8.2-1build3**. GRASS is a
  separately executed GPL-2.0-or-later program; no GRASS source is vendored.
- `vendor/rivgraph.py` extracts `algmap`, `dir_shortest_paths_links` and
  `guess_synthetic_slope` from the exact GRIT RivGraph fork revision
  `0a0a7135752894f451344bbf3521cc0f85f6f426`. Its BSD-3-Clause notice is retained
  in `vendor/LICENSE-RivGraph.txt`.
- `vendor/provenance.json` identifies normalized source hashes and modifications.
  `python -m river_flattening.vendor.refresh` downloads the pinned sources and
  reproduces the extract. Changes reject exact outlet-distance ties, preserve
  the shortest of parallel graph edges, and remove a module-namespace dependency.
- The continuity rule is an iterative adapter of RivGraph's rule: at a junction
  with one unknown link, all established links entering imply that the unknown
  link leaves (and vice versa). Our implementation checks conflicting evidence
  before propagating. It does not import RivGraph's recursive width assumptions,
  arbitrary fallback, or cycle-flipping routines.

This is a scoped integration of selected algorithms, not execution of the whole
GRIT or RivGraph workflow. It avoids installing their incompatible pinned Python
stacks. No GRIT/FABDEM data or GLO-30 elevations are used.

## Routing surface and nodata

`routing.py` creates a temporary analysis raster on the original geographic grid.
It interpolates holes at most 64 pixels from supported non-ocean CartoDEM data,
using GDAL's fill-nodata operation. Ocean remains NULL. Larger holes and missing
classification coverage remain NULL. Interpolation is processed in 512-pixel
blocks with a 72-pixel halo so block borders do not truncate the search.

This analysis raster is **never** an input to final elevation COGs. Final source
nodata stays `-255`. Neither routing interpolation nor its derived directions
creates valid output pixels. GRASS seeds drainage at internal NULL boundaries as
well as edges, so a hole must not be declared a sea mouth. Only GLO-30 class 1
supplies sea contacts.

A separate evidence raster excludes interpolated/missing pixels and their
eight-pixel neighbourhood from flow observations. Negative GRASS accumulation
marks incomplete boundary catchments and is excluded. Remaining drainage vectors
inside class 3 are assigned through the existing in-water ownership map, then
projected onto each local centreline tangent. The weight is
`min(log(1 + accumulation), 12)`. This is supporting evidence, not a measured
discharge or a calibrated confidence probability. Routing can still be wrong
where the source terrain or mask is wrong.

The initial runner deliberately limits a regional routing domain to 40 million
cells. Nationwide routing needs overlapping regional contexts and validation;
it must not silently run the entire national bounding rectangle in memory.

## Direction decisions

For each junction-to-junction or terminal reach:

1. Estimate a robust longitudinal bank trend from per-station CartoDEM bank
   quantiles with Theil–Sen slopes. Require at least three supported stations,
   90 m of evidence, an inferred drop of at least 3 m, and a 90% slope interval
   that excludes zero. At most 128 deterministically spaced stations are used
   for this direction calculation. The interval is a regression diagnostic,
   not a statement of absolute DEM accuracy.
2. Pass the resulting measured height trend into RivGraph's configurable
   attribute-slope guess. The synthetic outlet-distance default is not used
   as a substitute for CartoDEM height evidence.
3. Obtain GRASS votes along the reach. A station needs at least three eligible
   river samples and a projected mean vector magnitude of 0.15. A reach needs
   at least three voting stations and 75% agreement for a routing guess.
4. Generate RivGraph shortest-path guesses toward known ocean mouths, or a
   distinctly low terminal when supported terminal levels differ by at least
   3 m. Multiple sea mouths are allowed. Path geometry alone never establishes
   a direction.
5. Accept an ocean-mouth direction; otherwise accept a robust bank trend unless
   GRASS contradicts it, or accept agreement between GRASS and the outlet-path
   guess. Ocean/bank disagreements remain explicit diagnostics. A reach joining
   two ocean mouths stays unresolved as potentially tidal.
6. Propagate non-conflicting continuity. Invalidate directed cycles and links
   meeting at an unexplained interior source/sink. Do not randomly flip links.

These are conservative trial thresholds. Every reach receives bank, routing,
outlet-path, acceptance/rejection and profile-application diagnostics in
`directions.json`. The existing support gates (60% supported stations, 1,200 m
maximum bank-support gap) and station adjustment gates (95th percentile 15 m,
maximum 40 m) remain in force. Heights still come from the shared quantile profile
solver, not GRASS accumulation or a synthetic DEM.

## Partial-treatment boundaries

**2026-10-24 review: the 10 m gate is only a coarse failure diagnostic, not a
production acceptance standard.** Passing it does not make the remaining joins
acceptable. The v2 output remains an experimental pilot.

The first integrated pilot passed the native validity and downstream profile
checks but created 30 cardinal river-pixel pairs whose height step increased by
more than 10 m (maximum increase 16.80 m). This occurs where a fitted surface
meets retained original heights; monotonicity on accepted reaches alone cannot
detect it.

`boundaries.py` therefore checks neighbouring valid river pixels before COG
application, including contacts across sheet boundaries. If a final step exceeds
10 m and increases by more than 10 m relative to the coastal baseline, the whole
responsible shared profile group is retained unchanged. It repeats until no such
new steps remain, since a reverted group can expose another transition. Reverting
isolated water pixels would break the fitted cross-sections, so the gate operates
on whole groups. The 10 m threshold is a conservative trial acceptance limit,
not a claimed universal physical maximum river slope. Existing source jumps are
not repaired by this gate. Contact tables are stored per sheet for bounded-memory
scanning. The independent `audit.py` verifies the written COGs again.

The [USGS Lidar Base Specification 2025 rev. A](https://www.usgs.gov/ngp-standards-and-specifications/lidar-base-specification-digital-elevation-model-surface)
requires river surfaces level across the channel, descending with the surrounding
terrain and no higher than the adjacent banks. Its
[hydro-flattening reference](https://www.usgs.gov/ngp-standards-and-specifications/lidar-base-specification-appendix-2-hydro-flattening-reference)
describes elevation-constrained shoreline breaklines as the predominant method.
These are useful surface-shape requirements, not a direct accuracy specification
for 30 m CartoDEM. It does not give a universal permissible artificial join step.
Documented tidal/acquisition-level discontinuities are a separate exception.

For another production example, [WorldDEM editing, section 3.3](https://isprs-archives.copernicus.org/articles/XL-7-W3/1149/2015/)
describes downhill river levels in 50 cm increments, with operator input and
neighbouring-cell/visual QC. That increment is not a generic 50 cm seam tolerance.

The next method revision should impose compatible boundary levels in the fitted
surface, including transitions to retained reaches, and check the longitudinal
gradient. A shared physical boundary should have one modeled elevation to
numerical precision; adjacent raster cells may differ because of actual slope.
Reducing the 10 m rejection number alone does not supply this constraint.
Where bank/outlet/boundary evidence cannot support a coherent surface, retain the
larger affected reach and record the unresolved condition. Do not force arbitrary
ramps or infer that the source DEM's several-metre uncertainty permits seam jumps.
This stronger treatment has not yet been implemented; v2 retains the diagnostic
gate described above so its results remain reproducible.

## Reproduction and validation

Install the project's pinned Python requirements and a GRASS CLI. Preserve the
recorded GRASS version when reproducing a benchmark. From `egm2008/`, after
`build.py download` and `prepare` have verified/prepared the inputs:

```bash
python -m river_flattening.production run \
  --workdir work --label pilot-direction-v1 \
  --sheets cdnd43v.tif cdnd43w.tif
```

Individual stages are `geometry`, `routing`, `fit`, and `apply`. The explicit
`--bank-directions-only` option is an ablation test; it does not fabricate routing
evidence. Recipes record input, code and output hashes, versions and parameters.
An incompatible cached recipe fails rather than silently resuming.

```bash
python -m unittest -v river_flattening.test_profiles \
  river_flattening.test_geometry river_flattening.test_production \
  river_flattening.test_directions river_flattening.test_boundaries
```

The 19 river tests pass: numerical profiles, geometry, cross-sheet COG output,
multiple mouths, unequal braided branches, conflicting evidence, outlet-distance
ties, unresolved parallel reaches, routing-only nodata interpolation, boundary
contacts across sheets and whole-profile rejection at introduced steps. Real
pilot validation must additionally inspect coastline/bank errors, native pixel
preservation, sheet joins and unresolved branches. Successful direction inference
does not itself validate a published river surface.

References: [RivGraph](https://esurf.copernicus.org/articles/8/87/2020/),
[GRIT section 2.7](https://doi.org/10.1029/2024WR038308),
[GRASS r.watershed](https://grass.osgeo.org/grass-stable/manuals/r.watershed.html),
and the pinned source links in [REUSE_REVIEW.md](REUSE_REVIEW.md).
