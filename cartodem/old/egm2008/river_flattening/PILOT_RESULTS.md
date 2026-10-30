# D43V and D43W river pilot findings

This is the local record for reuse by other jobs. The method is described in the
[README](README.md); these results do not establish nationwide readiness.

**Current work uses WRIS river polygons.** The source-pin WRIS v2 experiment
treated 6.24% but imposed an unplanned exact raw-pixel constraint and is rejected
as the shared-boundary method. See [WRIS_INTEGRATION.md](WRIS_INTEGRATION.md).
The separate joint-fit pilot treats 22.14% and removes partial-treatment
interfaces; its remaining gradient concerns are recorded at the end of this file.
Earlier results below use
GLO-30 river geometry and have a different coverage denominator.

**The historical integrated direction v2 pilot reached 38.3% coverage.**
It uses GRASS terrain routing, selected RivGraph guesses and a boundary acceptance
gate. See the final section below. The historical two-sheet test superseded
the bounded pilot as evidence of the original method's coverage.
Its wider context exposed multiple ocean outlets and caused the current direction
rules to retain much more of the river unchanged. The 86.9% figure below applies
only to the original crop.

## Inputs and reproducible command

The pilot covers 75.08–75.24 E, 12.10–12.28 N near Nileshwaram/Valiyaparamba.
It uses existing EGM2008-corrected CartoDEM D43V elevations and the corresponding
GLO-30 WBM. Only the ocean and river classes affect the output. No GLO-30 elevation
values are read.

From `egm2008/`, after installing the dependencies, the local input paths used
for this pilot are:

```bash
python -m river_flattening.study \
  --dem corrected_cartodem/cdnd43v_egm2008.tif \
  --mask pilots/d43v/glo30-n12-e075-wbm.tif \
  --bounds 75.08 12.10 75.24 12.28 \
  --max-support-gap 1200 \
  --output river_flattening/pilots/d43v-v4

python -m river_flattening.plot_study river_flattening/pilots/d43v-v4
```

The output directory must be new. For repetition, choose a different output
directory or use the source snapshot in `pilots/d43v-v4/code/` as a Python package.
The original corrected TIFF uses NaN nodata; the study normalizes its invalid
pixels to `-255` and applies the ocean mask before fitting. This does not restore
any missing CartoDEM elevations. The production coastal baseline can also be
used as the elevation input.

`report.json` records input-window hashes, code hashes, dependency versions and
all trial parameters. Its `code/` directory contains the exact processing source
used. `stations.geojson` records bank support and fitted station heights.

## Result

| Check | Result |
| --- | ---: |
| Native window | 649 × 577 pixels |
| GLO-30 river pixels | 18,925 |
| River pixels with valid CartoDEM elevations | 17,381 |
| River pixels treated | 15,108 (86.9% of valid river pixels) |
| 95th percentile absolute elevation change | 10.81 m |
| Maximum absolute elevation change | 67.44 m |
| Land, class 2 and original nodata | Exactly unchanged |

The largest connected component uses 400 stations and 6,910 bank observations.
Its fitted profiles share confluence elevations and pass the downstream
monotonicity check. The solver converged in 20 iterations. Its 95th percentile
station adjustment relative to bank quantiles is 10.16 m; maximum 21.52 m.
These residuals describe agreement with noisy source observations, not accuracy
against an independent elevation reference.

The largest change is at **12.17638889 N, 75.13722222 E**: the corrected source
pixel is 72.3754 m and the candidate surface is 4.9385 m. The surrounding 5 × 5
window contains a pronounced spike, mostly 17–34 m around that centre, alongside
original nodata. This is consistent with the noisy river surface motivating the
work; it is not an independent measurement proving the fitted height. The
adjacent nodata pixels remain `-255`.

![EGM2008 baseline, candidate river surface and changes](pilots/d43v-v4/comparison.png)

## What the experiments established

The original bank samples are not a reliable standalone flow-direction source
here. Several reaches appear to rise toward the coast by 5–11 m. Selecting the
minimum observed bank pixel, or accepting every apparent local downhill trend,
would therefore be unreliable.

Within the bounded crop, the main component has a single terminal adjoining the GLO-30 ocean class.
The v4 experiment uses that geometry to establish the outlet, retains the bank
contradictions in the report, and fits a shared constrained surface. It does not
force the mouth to zero metres. All fitted heights come from CartoDEM evidence.

A gap in usable bank observations extends about 766 m along the network. The
default 600 m support limit correctly rejected this component. The v4 experiment
explicitly used 1,200 m to test interpolation through that gap; original missing
river pixels remain nodata. The later production driver defaults to this 1,200 m
trial limit; its suitability nationwide remains unvalidated.

Other components remain unchanged because there is insufficient direction
evidence, including a connected channel whose lowest terminal candidates differ
by only 2.23 m. Closed-loop fragments without a supported projection are also
retained unchanged.

An earlier OSQP trial did not converge for one small component. The current
kernel uses Clarabel, which passes all nine numerical/geometry tests and converges
for the supported components in this pilot. The code snapshot pins that choice.

## Remaining checks before the production river build

- Inspect the largest source-to-surface changes and curved-channel cross-sections.
- Validate representative inland rivers, broad/braided rivers and real sheet seams.
- Assemble complete connected networks before final per-sheet output; the pilot
  uses one bounded window and may truncate upstream context.
- Decide and record supported behavior for multiple mouths, tidal networks,
  ambiguous directions and gaps in bank observations.
- Run and validate the nationwide driver after resolving its direction policy.
  Final TerrainRGB and publication remain pending.

## Later two-sheet test: D43V and D43W

The `production.py` driver assembled both complete sheets before fitting, then
exported two pilot COGs. The retained reports are in [pilots/d43v-w/](pilots/d43v-w/).

| Sheet | Valid river pixels | Treated pixels | 95th percentile absolute change | Maximum absolute change |
| --- | ---: | ---: | ---: | ---: |
| D43V | 104,024 | 18,072 | 11.70 m | 38.96 m |
| D43W | 56,127 | 9,146 | 7.33 m | 19.02 m |
| Combined | 160,151 | 27,218 (17.0%) | — | — |

The network contains 54,903 skeleton nodes, 11,309 stations and 461 components;
174 components were fitted and 287 remained unchanged. The main coastal component
has multiple ocean outlets and was rejected by the prototype's single-outlet rule.
Other large components failed terminal-height separation, direction-tie or
bank-height consistency checks. These failures describe the current method's
limitations; they do not establish that the underlying rivers are unsolvable or
that every rejected component contains a source-data defect.

Both COGs use Float32 elevations, `NoData=-255`, and average elevation overviews
at factors 2, 4 and 8. Native values and validity were checked against the intended
river-stage output. These pilot files have not been published.

## Integrated directions: full D43V/D43W, 2026-10-24

**The 38.3% figure is experimental treatment coverage, not production-approved
coverage.** The 10 m transition gate only detects severe failures. Standard
hydro-flattening surface continuity needs stronger boundary constraints; see
[the subsequent standards review](DIRECTION_INTEGRATION.md#partial-treatment-boundaries).

[Direction integration](DIRECTION_INTEGRATION.md) documents the algorithms,
source revision, license, parameters and reproduction command. Local reports,
profiles and a small coastal comparison crop are retained in
[pilots/direction-v2/](pilots/direction-v2/). Large source and output rasters remain
in the build workspace. The label for reproduction is `pilot-direction-v2`.

| Method | D43V treated | D43W treated | Total / 160,151 valid river pixels |
| --- | ---: | ---: | ---: |
| Historical single-outlet method | 18,072 | 9,146 | 27,218 (17.0%) |
| Integrated direction v1 | 52,431 | 15,678 | 68,109 (42.5%) |
| Integrated direction v2, transition gate | 46,514 | 14,777 | 61,291 (38.3%) |

The network has the same 54,903 native skeleton nodes, 11,309 stations and
1,890 reaches. Directions are accepted for 432 reaches; 1,458 remain unresolved.
Accepted evidence includes robust bank slopes, GRASS/outlet-path agreement,
continuity and ocean contacts. Short reaches without sufficient observations
remain common. Twenty-nine reaches have conflicting bank/routing evidence;
thirteen are rejected for cycles and six for unexplained interior sources/sinks.

Of 227 connected groups of accepted reaches, v2 fits 206 and retains 21
unchanged: seven fail station-adjustment limits, four lack bank support and ten
fail the transition gate. There are 4,386 fitted stations. The historical main
coastal component with multiple sea mouths now has 12,029 of 26,950 valid pixels
changed (44.6%); it is partly resolved, not completely flattened.

V1 introduced 30 cardinal neighbour pairs whose final height step exceeded 10 m
and increased by more than 10 m. The maximum increase was 16.80 m. V2 retains
the ten responsible shared profile groups unchanged. Independent COG inspection
finds **zero** such pairs, with a maximum remaining increase of 9.01 m. This
trial acceptance test does not remove steps already present in the original DEM,
and a remaining sub-threshold transition is not automatically physically correct.

| Check | D43V | D43W |
| --- | ---: | ---: |
| 95th percentile absolute treated-pixel change | 10.00 m | 8.04 m |
| Maximum absolute change | 38.96 m | 20.17 m |
| Native land, lake and nodata pixels | Exactly unchanged | Exactly unchanged |
| Original validity mask | Exactly unchanged | Exactly unchanged |
| Stored downstream profile violations | 0 across the combined network | 0 |

Both COGs are Float32, `NoData=-255`, with average overviews at 2, 4 and 8.
GRASS 7.8.2 interpolated 6,032 pixels only on its disposable routing surface;
none became valid output pixels. The complete pinned-environment test suite
passes **33 tests**, including 19 river tests.

The scoped code integration and coastal regression are complete. Nationwide
routing still needs overlapping regional processing beyond the current
40-million-cell regional limit, plus representative inland/braided validation.
National river COGs, TerrainRGB packaging and publication remain unfinished.
The [noise-correction review](../DEM_NOISE_REVIEW.md) records a potential
upstream improvement; no denoising has been applied to these pilot elevations.

## WRIS joint-fit pilot, 24 October 2026

The selected river source is now WRIS. GLO-30 supplies ocean only. WRIS `01`,
blank and null features are River; `02` Sand and `03` Island are excluded.
The feature meanings were verified against the WRIS ArcGIS renderer, and the
full published Parquet was checksum-verified. See [WRIS_INTEGRATION.md](WRIS_INTEGRATION.md).

The first boundary revision (`pilot-wris-v2`) pinned each treated interface
pixel to its own original height within 1 mm. This was an implementation error
relative to the agreed shared-boundary approach. The 105 infeasible groups and
6.24% coverage reflect that rejected experiment; do not reuse its raw-height pins
as the intended method.

The replacement (`pilot-wris-v3`) fits connected sections together with shared
junction unknowns. Bank observations determine those levels. Unknown directions
can be proposed from an initial joint profile and are enforced and checked in a
second monotone fit. They remain explicitly marked model hypotheses. Incomplete
connected groups are retained instead of producing isolated fitted patches.

| Check | D43V | D43W | Total |
| --- | ---: | ---: | ---: |
| Valid WRIS river pixels | 182,803 | 298,505 | 481,308 |
| Treated river pixels | 69,535 | 37,010 | 106,545 (22.14%) |
| Original/ocean-masked nodata retained | 646,658 | 0 | 646,658 |
| Fitted-to-retained native river interfaces | 0 | 0 | 0, including sheet borders |
| Maximum absolute elevation change, m | 64.22 | 19.29 | 64.22 |
| 95th percentile absolute change, m | 11.63 | 6.70 | — |

All original nodata, nonriver and unsupported river samples are unchanged.
The final Float32 COGs have `NoData=-255`, lossless compression and averaged
numeric overviews. All stored profile edges are non-increasing downstream.
The complete current test suite passes **47 tests**, including the automatic
native COG audit. Source snapshots, reports, profiles and a small diagnostic
crop are retained in `pilots/wris-v3/`; the current auditor and its hash are also
recorded with the validation result.

There are 1,752 fitted connected groups out of 2,713 modelled groups. The applied
reaches comprise 538 with independent direction evidence and 1,534 whose
directions come from the joint-profile hypothesis. Rejections comprise 923 groups
with too little bank support, 30 needing excessive bank-level adjustment, one
with an unsupported section and seven with unresolved native connectivity.

The WRIS denominator differs from the historical GLO-30 denominator. WRIS includes
184,025 pixels that GLO-30 labels lake, mostly in D43W. GLO-30's lake class is
diagnostic only and does not override the user's selected WRIS river mask.
This is not a separate lake-flattening stage or a reason to compare the two
coverage percentages as if they described the same footprint.

**The surface is still experimental.** There are no added steps over the old
10 m diagnostic threshold, but that is not an acceptance criterion. Within
fitted groups, adjacent-pixel steps can increase by up to 9.59 m. At
12.85777778 N, 75.71083333 E, two neighbouring model positions change from
635.033/634.031 m to 624.434/635.029 m. Elsewhere, neighbouring pixels can be
assigned to skeleton positions hundreds of metres apart. These are within-model
gradient/projection concerns that shared junction heights alone do not solve.
Further slope and surface-projection validation is required before nationwide
publication; passing preservation and profile-order checks is insufficient.

The largest correction remains the known coastal-mask artefact at
12.17638889 N, 75.13722222 E: 72.375 m becomes 8.159 m. No new ground-truth water
elevation was established. The method and this result must not be described as
survey-validated.

COG SHA-256:

- D43V: `006efe809c5f73470b58854bc8d6a48b1bdcbe5cc1488fb1e37cfe52351c4476`
- D43W: `ece1c7c3b530eea01c9224c28bb5e96aadd6986516012c3e6edf97d324312d13`
