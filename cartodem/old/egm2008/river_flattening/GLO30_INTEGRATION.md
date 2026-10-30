# GLO-30 river selection across the CartoDEM footprint

Decision record: 24 October 2026. This supersedes WRIS as the production river
mask. Paths are relative to an arbitrary build workspace.

## Coverage and source choice

CartoDEM coverage extends beyond India. Valid prepared elevations were checked
at sample locations in Afghanistan, Pakistan, Nepal and Bhutan. All four
covering sheets also contain GLO-30 river-class pixels. These are point/sheet
checks, not a claim of complete coverage of every country.

| Sample location | Latitude, longitude | CartoDEM sheet | GLO-30 river pixels in the whole sheet |
| --- | --- | --- | ---: |
| Kabul, Afghanistan | 34.5553, 69.2075 | `cdni42j.tif` | 3,554 |
| Islamabad, Pakistan | 33.6844, 73.0479 | `cdni43n.tif` | 67,759 |
| Kathmandu, Nepal | 27.7172, 85.3240 | `cdng45b.tif` | 39,797 |
| Thimphu, Bhutan | 27.4728, 89.6393 | `cdng45f.tif` | 3,562 |

The India-focused WRIS selection cannot be assumed to supply river geometry
throughout these areas. Use one GLO-30 source throughout the DEM footprint,
including India, to avoid introducing a mask-source change at country borders.
There is no administrative boundary clipping or WRIS/GLO-30 blend.

GLO-30 is a global product, but the pinned public delivery's actual coverage must
still be checked. `build.py download` verifies classification availability:
missing GLO coverage is accepted only for an entirely nodata CartoDEM sheet.
Geometry checks every valid elevation against classification availability.
[Copernicus coverage](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM)
and the [public-delivery description](https://registry.opendata.aws/copernicus-dem/)
describe the product and distribution. Download URLs and checksums remain pinned
in `inputs.lock.json` and `reports/inputs.json`.

## Pixel selection and retained method

The [Copernicus DEM handbook, Table 8](https://dataspace.copernicus.eu/sites/default/files/media/files/2024-06/geo1988-copernicusdem-spe-002_producthandbook_i5.0.pdf)
defines WBM values as 0 no water, 1 ocean, 2 lake and 3 river.

- **Class 1:** clip valid ocean samples to `-255` during coastal preparation.
- **Class 3:** construct river geometry and apply supported fitted river heights.
- **Class 2:** preserve corrected CartoDEM lake heights; exclude lake pixels from
  river fitting and bank observations.
- **Class 0:** retain corrected heights; valid nearby pixels can supply bank
  observations. Land and islands remain outside river geometry.
- **Missing classification:** never infer land or river from absence.

Use native categorical masks with verified grid alignment and no overviews.
Do not dilate the selected river footprint or connect it across land. The
skeleton ownership and river projections remain confined to class 3.
GLO-30 supplies classifications only; no GLO-30 elevations are downloaded or
copied into the river fit. A 30 m classification cannot identify every narrow
channel or island; this source choice does not claim to remove that limitation.

The existing shared-junction solver is retained. Heights come from corrected
CartoDEM bank observations; connected reaches share fitted junction elevations,
and unsupported components stay at their corrected baseline. WhiteboxTools
provides direction evidence through the same staged breaching/D8 workflow.
Conditioning affects disposable routing rasters only. See
[NATIONAL_BUILD.md](NATIONAL_BUILD.md) and [WHITEBOX_TRIAL.md](WHITEBOX_TRIAL.md).

Original `-1024`, original `-255` and newly ocean-clipped samples all remain
`-255` nodata in COGs. River treatment never fills those holes. TerrainRGB nodata
remains the zero-height placeholder `(1, 134, 160, 254)`, with the same
zoom-dependent height quantization and lossless WebP encoding.

## Reproduction and cache separation

GLO-30 is the default `--river-source` in `run_build.py`, `build.py` and
`river_flattening.production`. WRIS remains an explicit comparison option.
The production download stage does not fetch WRIS unless requested.
The national COG tags, COG manifest, tile recipe, PMTiles attribution and generated
run notes identify the actual river source. Resuming after geometry with a
different source is rejected before detaching a worker.

From the source directory, using the pinned dependencies and precompiled
Whitebox installation:

```bash
python run_build.py --workdir /path/to/build --river-source glo30 \
  --label pilot-glo30-whitebox-v1 --sheets cdnd43v.tif cdnd43w.tif \
  --resume-from geometry --through cogs --detach
```

For a fresh full build omit the pilot label, sheets and stage limits. To reuse
verified coastal sheets after changing river source, archive the previous river
workspace and restart at geometry under an unused output label. The final
tiling path uses `final`; its previous river workspace must not be reused.

## Execution record

The first GLO-30 run below stopped at geometry after 185 sheets: constructing the
water-pixel graph for `cdnf46g.tif` exhausted memory. [ALLOCATION.md](ALLOCATION.md)
records the ready-made Whitebox replacement, successful full-worker test, and
new two-sheet pilot. Its predecessor's geometry checkpoints must be preserved
separately and rebuilt under the new allocation recipe.
That replacement build started at **12:12:29 UTC on 24 October 2026** under
`final`; the previous workspace is preserved as `glo30-scipy-final-20261024`.

The WRIS run was stopped after 14 routing contexts, before fitting or final COG
export. Its geometry, masks and direction samples are preserved under
`rivers/wris-final-20261024/`, with its log and status preserved separately.
The original inputs and all 543 prepared coastal sheets are reused. River
geometry and owner-based routing evidence must be rebuilt for the new footprint.

The full build restarted at **10:19:24 UTC on 24 October 2026** with
`run_build.py --river-source glo30 --resume-from geometry --detach`.
Its source snapshot is `reports/source-20261024T101924Z-glo30.tar.gz`, SHA-256
`664202f7d8d5ecc240aa3425e78bdef51a8581f3852449fdbceee5fd5a0e45ce`.
It verifies prepared inputs, rebuilds river geometry and routing, then continues
through fitting, COGs, numeric resampling, TerrainRGB and PMTiles. Live progress
is in `reports/build-final.json` and `logs/build-final.log`; a kickoff does not
establish completion of those stages.

The 57-test suite passed after the source switch, including unchanged lakes,
native nodata preservation, cross-sheet joint fitting, actual mask attribution
in a packaged PMTiles archive, and rejection of WRIS geometry during a default
GLO-30 restart.

The `pilot-glo30-whitebox-v1` production trial completed geometry, regional
Whitebox routing, shared component fitting and both COG exports in 127.14 seconds.
It treated **63,556 of 160,151 valid river pixels (39.7%)** and retained all
646,658 nodata samples. Both COGs declare `NoData=-255` and `RIVER_MASK=glo30`.
Native validity and every nonriver/nodata elevation were preserved exactly.
The audit found no partial-treatment interfaces or downstream station violations.
The maximum introduced adjacent-water jump was **6.5542 m**, so these remain
candidate surfaces requiring gradient review, not a claim of resolved joins at
every native pixel. Routing took 38.92 and 39.28 seconds for the two contexts.

Portable coverage checks, pilot measurements and source-change details are in
[pilots/glo30-whitebox-v1](pilots/glo30-whitebox-v1/). Full raster outputs remain
in `output/pilot-glo30-whitebox-v1/cogs/` in the build workspace.

The earlier 9.59 m introduced-step statistic and 22.14% treatment coverage belong
to the WRIS pilot. They must not be quoted as measurements of the new GLO-30
build. Earlier GLO-30 pilots also used older fitting/routing methods; their
percentages do not establish the current run's treatment coverage.
