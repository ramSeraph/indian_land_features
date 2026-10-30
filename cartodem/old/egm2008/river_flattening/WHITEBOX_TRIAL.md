# Whitebox routing and watershed trial

These original backend-validation records used WRIS geometry. The production
mask subsequently changed to GLO-30 class 3 across the full DEM footprint;
[GLO30_INTEGRATION.md](GLO30_INTEGRATION.md) records that change. The staged
Whitebox algorithms remain the same. The trial and production-equivalence
counts below describe the preserved WRIS run, not the new GLO-30 geometry.

Trial started 24 October 2026. This evaluates an alternative direction-evidence
backend; it does not change the agreed WRIS shared-junction height fit.

## Inputs and tools

The input is an EGM2008-corrected CartoDEM crop with GLO-30 ocean removed and
original nodata unified to `-255`. A sheet has 1024 pixels of neighbouring
context on each side, clipped to the prepared mosaic bounds. The trial can
prepare this crop directly from the corrected sheets and classification mosaic.
**No GRASS installation, process or output is required.** An existing GRASS
result can optionally supply a comparison baseline on the exact same DEM.

Small holes may be interpolated on the disposable routing DEM using the existing
`analysis_window` preparation. An independent evidence mask excludes original
holes and their uncertainty buffer. No provisional interpolation or subsequent
Whitebox conditioning is transferred to published COGs.

Whitebox 2.4's GeoTIFF reader does not support floating-point predictor 3, used
by our compressed analysis TIFF. The runner writes a predictor-1 copy and checks
every elevation sample, the grid and nodata for exact preservation before routing.
This changes only the TIFF encoding.

The successful workflow runs `BreachDepressions`, `D8Pointer`, and
`D8FlowAccumulation --pntr` in separate processes. Accumulation is in **cells**,
without logarithms or display clipping. Breaching uses upstream defaults and
writes its routing surface as Float64. The standard D8 pointers follow steepest
descent on that surface.

The initial [`FlowAccumulationFullWorkflow`](https://github.com/jblindsay/whitebox-tools/blob/28ed01f506f2024cfb52db85348b887f2ded4e5d/whitebox-tools-app/src/tools/hydro_analysis/flow_accum_full_workflow.rs)
trial used an aspect-aligned D8 method. It exhausted memory on the same crop.
It remains available through `--workflow full` for reproduction. The staged
method is a different algorithm, not a claim of identical results obtained by
splitting that combined tool.

[`Watershed`](https://github.com/jblindsay/whitebox-tools/blob/28ed01f506f2024cfb52db85348b887f2ded4e5d/whitebox-tools-app/src/tools/hydro_analysis/watershed.rs)
then delineates the upstream area of a selected reliable WRIS river pixel using
those pointers. Its cell count must match the accumulation at that point, allowing
only the accumulation raster's Float32 representation. This checks consistency,
not agreement with a surveyed catchment boundary.

## Missing upstream terrain

Whitebox's pointer codes differ from the internal compass encoding. The adapter
maps all eight explicitly and rejects unknown codes. Valid cells bordering nodata
or the crop edge seed an uncertainty mask, which follows the **actual Whitebox
pointers downstream**. Those cells cannot supply reliable river-direction samples.
The rule is conservative at ocean edges. The existing source-validity and eight-
pixel context buffers also apply.

Whitebox includes an MIT-licensed
[`EdgeContamination` source tool](https://github.com/jblindsay/whitebox-tools/blob/28ed01f506f2024cfb52db85348b887f2ded4e5d/whitebox-plugins/src/edge_contamination/main.rs),
but its D8 mode derives steepest-descent directions from a conditioned DEM. That
is not necessarily the aspect-aligned pointer field of the initial combined
workflow. The adapter traces supplied pointers for either workflow. The resulting exclusions are not claimed
to reproduce GRASS's negative-accumulation classification exactly.

## Precompiled installation and reproduction

The default installer is [install_whitebox_binary.py](../install_whitebox_binary.py).
It downloads the precompiled Linux ARM64 package from
[conda-forge](https://anaconda.org/conda-forge/whitebox_tools/files?channel=main):
`whitebox_tools-2.4.0-py314hfe60d44_2.conda`, SHA-256
`e3e680dc957a87fa99ed73f9cd0fc9028866557684af28b09ac9537e76dc2102`.
It also supplies a checksum-pinned private libgcc 14.3 runtime, retains both
packages' licenses, verifies the executable, and records provenance. It uses the
standalone CLI without installing Conda or the package's Python wrapper.
The Python 3.14 project environment supplies Zstandard decompression.

From the `egm2008` source directory, after the normal preparation and WRIS geometry
stages:

```bash
python install_whitebox_binary.py --workdir /path/to/build
python -m river_flattening.whitebox_trial \
  --workdir /path/to/build --sheet cdnc43d.tif \
  --geometry-label final --label whitebox-prebuilt-v1
```

Use the pinned project Python environment. `--context /path/to/crop` optionally
reuses `analysis.tif` and `evidence.tif` instead of preparing them. Separately,
`--grass-reference /path/to/reference` enables comparison with existing
`drainage.tif` and `accumulation.tif`; it never runs GRASS. When the reference DEM
is present its SHA-256 must equal the trial input's SHA-256.

The runner records execution time, observed peak memory, inputs and output hashes,
native grid and nodata checks, river evidence counts, and watershed consistency
in `rivers/LABEL/report.json`. Its current RSS monitor samples every 20 ms and
terminates a Whitebox process observed above 1100 MiB. This sampling is not a
kernel-enforced memory bound. Failed or incomplete runs are not accepted. Each
rerun uses a fresh label. Results do not alter the production routing cache.

## Source-build fallback

The bootstrap [install_whitebox.py](../install_whitebox.py) pins WhiteboxTools
Open Core 2.4.0 at commit `28ed01f506f2024cfb52db85348b887f2ded4e5d`, verifies the
source archive and upstream Cargo lockfile, and builds with Rust 1.79.0 and
`cargo build --release --locked --jobs 1`. This older compiler is compatible with
upstream's locked `time` dependency. All downloads, toolchains and build outputs
stay in the chosen workspace. The checksum-pinned bootstrap currently supports
Linux aarch64 and requires a working system C linker.

The original workspace lockfile failed Cargo validation for `tsp-rs 0.1.0`, an
unused travelling-salesman plugin dependency. Its published registry checksum
differs from the upstream lock. The build excludes the plugin and GUI workspaces,
then uses [whitebox-core.Cargo.lock](../whitebox-core.Cargo.lock), containing only
the 141 packages retained by Cargo for the core workspace. Every retained version
and checksum is unchanged from upstream. No replacement `tsp-rs` code is accepted and
no hydrology source is edited. Both original and focused lock hashes are recorded.

Compiling the complete command-line tool registry exhausted memory on a 2 GiB
system: the kernel killed `rustc` before DEM processing began. The bootstrap
therefore registers eight hydrology tools: the full flow workflow, Watershed,
Basins, D8Pointer, D8FlowAccumulation, BreachDepressions,
BreachDepressionsLeastCost and FillDepressionsWangAndLiu. Only the two registry
modules change. Each selected algorithm file is verified against the pinned
source archive; registry and algorithm hashes are recorded. The smaller build
completed in 24 seconds using the already-built dependencies. The first compiler
RSS sampler missed children of Cargo worker threads; it has been corrected and
the original 36.58 MiB compiler observation must not be treated as its peak.

The fallback trial selection is `--tool-installation tools/whitebox-2.4.0` after
running `python install_whitebox.py --workdir /path/to/build`. Precompiled
installation is preferred; compiling is unnecessary for the supported target.

## Execution record

The four adapter checks passed: compass directions, unknown-code rejection,
downstream propagation from missing upstream terrain, and termination on a cycle.
The staged real-sheet test on `cdnc43d.tif` completed on the exact same analysis
raster as the saved reference. Its crop contains 31,899,904 cells, of which
14,305,645 are valid on the disposable routing surface.

| Tool | Seconds | Observed peak RSS, MiB |
| --- | ---: | ---: |
| BreachDepressions | 20.10 | 801.03 |
| D8Pointer | 4.59 | 489.89 |
| D8FlowAccumulation | 6.75 | 494.07 |
| Watershed | 2.94 | 1007.81 |

Whitebox supplied accepted direction samples at 122,860 of 130,258 WRIS river
pixels, versus 123,718 for GRASS under their respective uncertainty exclusions.
Of 122,037 jointly accepted samples, 89,922 have identical directions and 20,158
differ by 45 degrees. These are agreement and coverage measurements, not accuracy
against independently observed flow directions.

The selected pour point at 11.51611111 N, 75.64305556 E has 939,907 upstream cells;
the number of cells labelled by `Watershed` exactly matches accumulation.
Grid and nodata checks passed. Temporary breaching lowered 1,225,836 analysis
samples by up to 58.0651 m; these are routing-only changes. No published elevation
sheet or river-profile fit was changed by this trial.

The precompiled conda-forge binary repeated the entire trial successfully.
Its six output TIFFs have **identical SHA-256 hashes** to the source build:
conditioned DEM, pointer, accumulation, uncertainty, pour point and watershed.
All comparison statistics are also identical. Its measured times were 17.53 s,
4.79 s, 6.62 s and 3.40 s respectively for the four tools; observed peak RSS was
804.52, 496.84, 509.77 and 1011.33 MiB. Compiling is unnecessary on this target.
Portable reports and the equivalence record are in
[pilots/whitebox-routing](pilots/whitebox-routing/).

The production adapter in `whitebox.py` uses the same three staged routing
commands and uncertainty mask. `regional.py` converts pointer codes to the
existing compass convention and samples only the owning sheet's river pixels.
The joint shared-junction fit is unchanged. The full 55-test suite passed after
this integration. The national runner defaults to this precompiled backend;
GRASS remains optional for historical comparisons. Routing caches include the
backend, binary/runtime/settings and adapter-code hashes.

The resumed national run's first context also matched the independent trial
exactly: 122,860 accepted samples across 31,358 owners, with identical vectors,
weights, counts and accumulation TIFF checksum. Reproduce that comparison with
`python -m river_flattening.validate_whitebox_routing --workdir /path/to/build`.
Its [portable verification record](pilots/whitebox-routing/production-equivalence.json)
and the [restart record](pilots/whitebox-routing/restart.json) preserve the result.

The initial combined workflow was killed by the kernel for memory exhaustion
after writing its conditioned DEM. Its 200 ms RSS sampler last observed
1385.47 MiB, below its then-1500 MiB cutoff; this did not prevent the subsequent
allocation burst. The successful staged run used the lower, faster monitor.

The existing national GRASS stage stopped on its first routing region. The saved
DEM subsequently ran successfully through the same worker, including a repeat
with the original temporary-directory setting. Therefore that setting has **not**
been established as the cause. The Whitebox replacement proceeds independently
of that unresolved failure.
