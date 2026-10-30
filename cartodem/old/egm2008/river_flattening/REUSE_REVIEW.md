# Existing river-direction code: reuse review

Reviewed 2026-10-23. This file records the source-code review that preceded the
integration. The subsequent [implementation record](DIRECTION_INTEGRATION.md)
identifies the selected RivGraph functions actually reused, GRASS version,
adapters and tests. The latest integrated two-sheet run treats 38.3% of valid river
pixels, versus 17.0% with the historical method, after a transition acceptance
gate. The coastal regression passes; nationwide validation remains unfinished.

## Recommendation

Evaluate **GRASS `r.watershed` for terrain-derived direction/accumulation evidence**
and **RivGraph's network-direction machinery for branching channels and multiple
outlets**. GRIT provides a published example of combining these approaches.
Reuse those components through an adapter to our network and profile solver.
The full GRIT data-production workflow is not required for this cartographic DEM.

Our inputs would remain EGM2008-corrected CartoDEM and GLO-30 ocean/river masks.
Direction inference must produce evidence and diagnostics for the existing
height-fitting stage. It does not authorize changing land, lake class 2 or the
original nodata footprint. GLO-30 elevations are not required.

## Sources and exact revisions inspected

| Component | Revision | License | Relevant source |
| --- | --- | --- | --- |
| RivGraph upstream | `f8745a645d822377c75e16f97ae995da32587695` | BSD-3-Clause | [directionality.py](https://github.com/VeinsOfTheEarth/RivGraph/blob/f8745a645d822377c75e16f97ae995da32587695/rivgraph/directionality.py), [delta directionality](https://github.com/VeinsOfTheEarth/RivGraph/blob/f8745a645d822377c75e16f97ae995da32587695/rivgraph/deltas/delta_directionality.py), [license](https://github.com/VeinsOfTheEarth/RivGraph/blob/f8745a645d822377c75e16f97ae995da32587695/LICENSE.txt) |
| GRIT orchestration | `1e9980bc3b8f62112a9754346c0b086e7154a9b4` | Apache-2.0 | [multi_network_rivgraph.py](https://github.com/mwort/grit/blob/1e9980bc3b8f62112a9754346c0b086e7154a9b4/src/multi_network_rivgraph.py), [river_network_utils.py](https://github.com/mwort/grit/blob/1e9980bc3b8f62112a9754346c0b086e7154a9b4/src/river_network_utils.py), [license](https://github.com/mwort/grit/blob/1e9980bc3b8f62112a9754346c0b086e7154a9b4/LICENSE) |
| GRIT's RivGraph fork/submodule | `0a0a7135752894f451344bbf3521cc0f85f6f426` | BSD-3-Clause | [RiverNetwork](https://github.com/mwort/RivGraph/blob/0a0a7135752894f451344bbf3521cc0f85f6f426/rivgraph/networks/__init__.py), [network_directionality.py](https://github.com/mwort/RivGraph/blob/0a0a7135752894f451344bbf3521cc0f85f6f426/rivgraph/networks/network_directionality.py), [license](https://github.com/mwort/RivGraph/blob/0a0a7135752894f451344bbf3521cc0f85f6f426/LICENSE.txt) |
| GRASS | main resolved to `51551cd8568669db60b666ef25d0afcd6e622367` during review | GPL-2.0-or-later, subject to component notices | [r.watershed source](https://github.com/OSGeo/grass/tree/51551cd8568669db60b666ef25d0afcd6e622367/raster/r.watershed), [COPYING](https://github.com/OSGeo/grass/blob/51551cd8568669db60b666ef25d0afcd6e622367/COPYING) |

GRIT's published data have a separate
[CC BY-NC 4.0 license](https://github.com/mwort/grit/blob/1e9980bc3b8f62112a9754346c0b086e7154a9b4/docs/GRIT/LICENSE.md).
That is distinct from the software licenses above. No GRIT/FABDEM dataset was
downloaded or adopted during this review. If code is incorporated, retain its
license and attribution and identify our modifications.

## What is actually reusable

RivGraph accepts explicit inlet/outlet lists. Its direction routines combine
network connectivity, shortest paths, channel geometry and already-established
directions. Useful functions include `set_inletoutlet`, `dir_shortest_paths_nodes`,
`dir_shortest_paths_links`, `set_continuity`, `set_by_known_flow_directions` and
cycle checking/repair. These cover topology that our single-outlet prototype
currently rejects.

The GRIT fork adds `RiverNetwork`, `set_width_continuity`,
`set_by_guess_agreement`, `drainage_area_gradient` and `guess_synthetic_slope`.
The class consumes link/node tables, which makes it a candidate for adapting
our already-extracted network rather than repeating mask extraction.
`multi_network_rivgraph.py` demonstrates processing connected components and
exporting link-direction and continuity diagnostics.

The [GRIT paper, section 2.7](https://doi.org/10.1029/2024WR038308)
describes combining accumulation, elevation, outlet distance, channel angles,
widths and network continuity. It also reports unresolved cycles requiring
removal or manual direction assignment. It is evidence for the approach, not
a guarantee that every network can be resolved automatically.

## Differences between the description and the inspected code

The GRIT input builder `rivgraph_links_input` computes `elevation_start` and
`elevation_end` as medians over the first and last halves of each sampled reach,
then stores their difference as `elevation_gradient` (lines 666–670).
These are sampled along the line, not specifically from surrounding banks.

However, the pinned fork's `set_inexact` does **not** directly consume those
elevation fields. It uses drainage-area gradients, outlet-distance gradients,
main-channel drainage-area gradients, angular agreement, widths and shortest
paths. `guess_synthetic_slope` accepts a configurable node attribute, but its
default `synthetic_elevation` is outlet distance in the GRIT input builder.
It must not be described as already using our CartoDEM bank elevations.
Adding robust CartoDEM reach-height evidence requires an explicit adapter or
documented modification to the direction-guess stage.

The upstream delta routine also uses a **synthetic** DEM made from inlet/outlet
locations. This supplies a geometric direction guess, not measured elevation.

The fork marks remaining unset directions with algorithm code 29 (`random`)
and marks them certain before cycle repair. At that point it retains their
existing orientation; this is not measured evidence. A CartoDEM integration must
retain unresolved status instead of treating that flag as confidence. Other
routines also contain random sampling, so reproducibility needs explicit handling
if those paths are used.

Use the exact GRIT submodule revision above. The fork's default `master` branch
does not contain `rivgraph/networks/network_directionality.py`; fetching it there
returned 404. Installing ordinary upstream RivGraph is therefore insufficient
to obtain the GRIT extension.

## GRASS role and limitations

[`r.watershed`](https://grass.osgeo.org/grass-stable/manuals/r.watershed.html)
uses least-cost terrain routing, offers multiple-flow accumulation and a disk
mode with a memory limit. Its `drainage` raster still records only the strongest
direction per cell. It can supply supporting terrain evidence; that raster alone
does not preserve all branches of the mapped river network.

Its initialization seeds raster edges and cells beside missing data. Original
CartoDEM holes must therefore be distinguished from real sea outlets when
preparing an analysis surface. Any provisional interpolation for routing would
belong only to a separate analysis surface, with provenance; final COG nodata
must remain unchanged. Cropped boundaries also require care. See
[init_vars.c](https://github.com/OSGeo/grass/blob/51551cd8568669db60b666ef25d0afcd6e622367/raster/r.watershed/ram/init_vars.c).

## Integration checklist from the source review

1. Convert our connected reaches, metric lengths/widths and verified ocean
   contacts to the expected tables. Do not assume every non-ocean endpoint is
   an inlet: mask gaps, lake contacts and crop boundaries can create endpoints.
2. Choose whether the first pilot uses regional GRASS direction/accumulation
   evidence, a CartoDEM height-trend adapter, or both. The existing GRIT recipe
   requires accumulation fields we do not currently compute; do not fabricate
   them just to satisfy the interface.
3. Preserve direction provenance, unresolved links and cycle/continuity checks.
   Their heuristic `certain` flag is not a calibrated probability.
4. Pass supported directions into our shared profile fit and retain native
   validity during application. Do not transfer any routing-only elevation
   alterations into the final COGs.
5. Test the complete coastal network that defeated the two-sheet prototype,
   followed by inland and braided examples. Report coverage, changes and failed
   checks; successful execution alone does not establish correct directions.

Dependency compatibility also needs a separate environment or a scoped port.
The fork declares old NetworkX/Shapely constraints; the inspected current upstream
package pins NumPy below 2 and older SciPy/Rasterio than our build. No package
installation, direction-code replacement or production rerun was performed during
the source-review phase itself. The later integration vendors only selected,
attributed functions rather than installing the incompatible complete stack.
