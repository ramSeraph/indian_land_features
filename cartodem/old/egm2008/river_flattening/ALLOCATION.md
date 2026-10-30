# River pixel allocation without the large water graph

Method and benchmark record: 24 October 2026. Paths are relative to an arbitrary
build workspace. This changes water-to-centreline assignment; the agreed shared
junction fit, direction evidence, datum conversion and nodata rules remain.

## Failure and existing-code choice

The GLO-30 build stopped after 185 of 543 geometry sheets, before routing. The
isolated worker for `cdnf46g.tif` was killed by the kernel for exhausting memory.
Its 4,624 × 4,624 context contains 2,749,455 river pixels and 156,133 retained
centreline pixels. The previous water graph had 21,128,122 directed neighbour
entries. Constructing temporary coordinate arrays before conversion to CSR
caused the allocation failure; process isolation alone could not fix it.

Use the already installed, checksum-pinned **WhiteboxTools 2.4.0 `CostDistance`
and `CostAllocation`** implementations. Their
[official documentation](https://jblindsay.github.io/wbt_book/available_tools/gis_analysis_distance_tools.html#costdistance)
describes raster cost propagation and allocation from source labels and
back-links. This removes the explicit graph of every water pixel; no new
shortest-path implementation is required. The smaller skeleton graphs still
use the existing SciPy code.

Other existing code was considered:

- [scikit-image `MCP_Geometric`](https://scikit-image.org/docs/stable/api/skimage.graph.html#skimage.graph.MCP_Geometric)
  supplies implicit raster paths and anisotropic sampling, but its public result
  needs an additional traceback-to-owner conversion. It was not benchmarked.
- scikit-image's `pixel_graph` also constructs explicit neighbour arrays, so it
  does not remove the problematic representation.
- [GDAL Polygonize](https://gdal.org/en/stable/programs/gdal_polygonize.html)
  produces footprint polygons. Those boundaries do not supply the centreline
  or the required assignment of water pixels to it. Vector polygons are not
  required for this raster-mask workflow.

## Adapter and preserved geometry

`allocation.py` writes two disposable rasters: cost 1 inside river class 3, nodata
outside, and one unique positive Int32 label per retained skeleton pixel. The
tools perform eight-neighbour least-cost allocation. Paths can use diagonal
river neighbours, as in the original graph, but cannot traverse land cells.
Separate channels and island interiors remain barriers.

Temporary raster cell sizes equal the existing window's `dx` and `dy` in metres,
including the cosine-of-latitude correction. Their EPSG:3857 tag is only a
disposable metric plane, not a reprojection or a claim about geographic position.
No published raster passes through that plane. Distances and diagonal lengths
therefore implement the existing graph metric rather than treating longitude
and latitude degrees as equal ground distances.

`CostAllocation` labels are converted back to native skeleton pixel IDs and then
to the common mosaic's global pixel IDs. The 512-pixel halo, spur pruning, bank
selection, shared-junction solver and projection remain. Every output block is
checked for grid identity, exact river-mask coverage, valid source labels and
unchanged seed identity. An unassigned river pixel fails the stage.

Geometry recipes record the backend, adapter hashes, binary/runtime hashes and
Whitebox settings. Per-sheet input/output hashes, tool timings and logs remain
under `rivers/<label>/allocation/<sheet>/`. Disposable TIFFs are removed only
after successful validation; a failed operation retains its scratch files.
The tools have the existing 1,100 MiB sampled RSS guard, which is not a hard
kernel memory limit.

The pure `extract_window` helper retains SciPy for small reference tests;
production geometry defaults to Whitebox. `--allocation-backend scipy` on the
geometry command is an explicit comparison option, unsuitable for this large
failing context. GRASS is not required by the new allocator or current national
routing path.

## Validation and measured limits

| Check on the failing context | Time | Peak observed RSS |
| --- | ---: | ---: |
| `CostDistance` alone | 3.35 s | 691.75 MiB |
| `CostAllocation` alone | 1.89 s | 644.33 MiB |
| Complete isolated geometry worker | 25.83 s | 1,145.03 MiB, process tree |

The first two rows exclude skeleton extraction and raster preparation. The last
row includes the actual worker's raster reads, skeleton pruning, allocation,
bank sampling and geometry export on the same full mosaic context. It completed
successfully. Sampling intervals were 20 ms for individual tools and 50 ms for
the worker tree; short-lived peaks can be missed. Small test processes that
finished between samples have reported RSS zero, meaning **unmeasured**, not
zero memory consumption.

Four synthetic comparisons cover an island, separated parallel channels,
equal-distance ties and a diagonal-only connection, including unequal cell
dimensions. Assigned sources have minimum path cost according to the old SciPy
graph. Float32 distance-file rounding was below 0.00013 m. The production
adapter also passed a two-window halo/global-owner check. All 61 unit tests
passed, including nodata preservation, cross-sheet fitting and transport quoting.

The complete `pilot-glo30-allocation-v1` two-sheet production pilot took 132.74 s.
It treated 63,556 of 160,151 valid river pixels and preserved all 646,658 nodata
pixels and all nonriver elevations exactly. Against the preceding
`pilot-glo30-whitebox-v1`:

- Water owners differed at 107 pixels; bank-owner assignments differed at 53.
- `cdnd43v.tif` has identical elevation samples.
- `cdnd43w.tif` differs at 469 river pixels, by at most 0.2609863 m.
- Treatment counts, partial-treatment-interface counts and downstream station
  violations are unchanged. The maximum introduced adjacent-water jump remains
  **6.5542 m**, so the surfaces remain candidates for gradient review.

Equal-cost source selection is not guaranteed identical between SciPy and
Whitebox. The synthetic tests demonstrate that difference; the real-data pilot
measures its practical consequences. This is not a byte-identical substitution,
and unchanged pilot coverage does not establish nationwide correctness.

Portable measurements are saved locally in [pilots/allocation-v1](pilots/allocation-v1/).
The two source rasters, tool outputs and detailed logs for the full allocation
trial remain under `rivers/allocation-f46g-v1/` in the build workspace.

## Reproduction

Use the pinned requirements and `install_whitebox_binary.py`. From the source
directory, against an existing verified geometry context:

```bash
python -m river_flattening.allocation_trial synthetic \
  --workdir /path/to/build --label allocation-synthetic-v2
python -m river_flattening.allocation_trial prepare \
  --workdir /path/to/build --label allocation-f46g-v1 \
  --geometry-label final --sheet cdnf46g.tif
python -m river_flattening.allocation_trial run \
  --workdir /path/to/build --label allocation-f46g-v1
python -m river_flattening.allocation_trial worker \
  --workdir /path/to/build --label allocation-worker-f46g-v1 \
  --geometry-label final --sheet cdnf46g.tif
python run_build.py --workdir /path/to/build --river-source glo30 \
  --label pilot-glo30-allocation-v1 --sheets cdnd43v.tif cdnd43w.tif \
  --resume-from geometry --through cogs --detach
# After the pilot completes, compare against its saved preceding pilot:
python -m river_flattening.allocation_comparison --workdir /path/to/build
```

Use unused trial labels. After archiving the old `final` context, pass its archive
label as `--geometry-label`. Trials do not modify production checkpoints.

The 185 SciPy geometry checkpoints must not be mixed with the new backend.
Preserve the failed river workspace, log and status under separate archive names,
then restart at geometry with `run_build.py --workdir /path/to/build
--river-source glo30 --resume-from geometry --detach`. Verified inputs and all
543 prepared coastal sheets are reused. Geometry, owner-based routing evidence
and downstream products are rebuilt with the new recipe.

## Restart record

The failed 185-checkpoint run is preserved as
`rivers/glo30-scipy-final-20261024/`, with correspondingly named archived log and
status files. The new full build started at **12:12:29 UTC on 24 October 2026**
under `final`, beginning at geometry. This records a restart, not completed
COGs or tiles. Current progress remains in `reports/build-final.json` and
`logs/build-final.log`.

The frozen source archive is `reports/source-20261024T121223Z-allocation.tar.gz`,
SHA-256 `abf0984fa27ed4ac971783313be758ce7a77bce2b746794f80f8a449198575ba`
(303 source/documentation/record files). Private connection configuration is
excluded. [The portable switch record](pilots/allocation-v1/allocation-backend-switch.json)
contains the archive paths, validation reference and preserved-input counts.
