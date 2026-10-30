# CartoDEM EGM2008 TerrainRGB build plan

Rebuild the CartoDEM terrain tiles from the original V3R1 sheets, correcting elevations to EGM2008 and applying the GLO-30 coastal mask before tiling. Preserve missing elevations throughout processing, store a uniform nodata value in the intermediate sheets, and encode nodata as exactly zero metres with alpha 255 in lossless WebP display tiles. Retain the previous zoom-dependent height quantization. COGs retain declared NoData=-255; display tiles do not preserve a distinguishable nodata marker. For existing completed COGs, use [the tiling-only workflow](TILING_SETUP.md).

The final product includes a river-treatment stage, followed by publishable native-grid COGs and TerrainRGB. River components without sufficient evidence remain unchanged and are reported. The [local river method record](river_flattening/README.md) contains the implementation, references, pilot results and remaining validation work.

## Reproducible workspace

Choose a working directory with enough space for the source sheets, intermediate rasters and packaged tiles. Scripts take this directory as an argument; downloads, conversion and packaging must run without host-specific paths or connection details. Use this directory layout:

```text
WORKDIR/
    code/                       Conversion and tiling scripts
    inputs/cartodem/             Original CartoDEM sheets
    inputs/glo30-wbm/            GLO-30 ocean and river classification
    inputs/wris/                Optional historical WRIS comparison inputs
    inputs/geoid/               Pinned EGM2008 geoid grid
    sheets/coastal/              Corrected and coast-masked sheets
    rivers/final/               Shared river geometry, fitted profiles and reports
    mosaics/                    VRTs referencing prepared sheets
    tiles/                      TerrainRGB tile output
    output/final/cogs/           Final publishable elevation COGs
    output/final/                Packaged PMTiles
    reports/                    Input manifest and validation results
    logs/                       Processing logs and completion records
    tmp/                        Temporary files on the external volume
```

At execution time, check free space, available memory, CPU count and software versions, then choose bounded worker counts. Keep temporary data under the working directory. The download stage must fetch missing inputs and verify existing files against the recorded inventory and checksums. Record dependencies and commands so the same build can be reproduced in another environment.

## Inputs

| Input | Purpose |
| --- | --- |
| [Original CartoDEM 30 m V3R1 release](https://github.com/ramSeraph/indian_land_features/releases/tag/cartodem-30m-v3r1-orig) | Source elevations and sheet index |
| `us_nga_egm08_25.tif` | Spatially varying EGM2008 geoid undulation |
| GLO-30 water-body masks covering the source sheets | Ocean class 1 and river class 3; lake class 2 excluded from river treatment |

Use the geoid grid and checksum already specified in the [README](README.md#the-geoid-grid-and-correction-sign). Record input URLs, versions, checksums and processing options in a build manifest. Use relative paths in manifests and generated VRTs. GLO-30 elevation rasters are not inputs to the water treatment; water heights come from CartoDEM. Use GLO-30 ocean and river classes consistently across the whole DEM coverage, which extends outside India. Heights outside the river footprint retain CartoDEM's existing treatment, apart from datum correction and ocean clipping. There is no additional lake-flattening stage.

## Processing sequence

```text
Original CartoDEM sheets
    -> identify source validity
    -> correct valid elevations to EGM2008
    -> apply GLO-30 ocean mask
    -> write Float32 sheets with NoData=-255
    -> fit shared river profiles across sheets using CartoDEM-derived heights
    -> apply supported profiles and export per-sheet Float32 COGs
    -> mosaic final COGs through a VRT
    -> resample elevations and validity for each zoom
    -> quantize valid elevations and encode TerrainRGB
    -> write lossless RGBA WebP and package PMTiles
```

### Identify nodata before correction

Treat both original values, **`-1024` and `-255`**, as nodata even when the source TIFF declares neither. Also honor existing raster masks and nonfinite samples.

Keep validity separate from numeric elevations during processing. Zero and other negative source heights remain valid unless explicitly masked. Do not subtract the geoid from a source sentinel or infer validity by comparing corrected elevations with the original sentinel values.

### Correct elevations to EGM2008

Use the existing source-height assumption: original CartoDEM samples are WGS84 ellipsoidal heights in metres. Bilinearly sample the pinned geoid grid at the original pixel centres, then calculate:

```text
H = h - N
```

Here `h` is the original elevation, `N` is geoid undulation, and `H` is the resulting EGM2008 elevation. Preserve the source CRS, dimensions, affine transform and pixel alignment. Fail the sheet if the geoid is unavailable at an otherwise-valid source pixel.

Apply this correction once. Later Web Mercator reprojection must change horizontal coordinates without applying another vertical transformation.

### Apply the coastal mask

Use the GLO-30 water-body mask's ocean class, value **1**, to mark ocean pixels nodata. Verify the mask's georeferencing against the source sheet; align categorical masks using nearest-neighbour sampling when grids differ.

The resulting validity is:

```text
prepared_valid = source_valid AND NOT glo30_ocean
```

Require mask coverage for the area being processed. Missing mask data must not silently become land or ocean. Retain valid CartoDEM elevations outside the ocean class, including lakes and rivers during the initial build.

Where GLO-30 identifies land but CartoDEM already contains nodata, retain nodata. A better coastline mask cannot recover elevations removed by the original CartoDEM mask. Report these locations separately from valid source pixels newly removed by coastal masking.

### Write intermediate sheets

Write tiled, compressed **Float32 GeoTIFFs**, using DEFLATE with floating-point predictor 3. Preserve native sheet geometry and declare **`NoData=-255`**. Float32 retains the fractional heights introduced by geoid correction.

Both original sentinel values and newly masked ocean pixels become `-255`. Retain the original per-sentinel counts in reports, since the stored output no longer distinguishes them. Record the EGM2008 datum, source and mask provenance, and processing settings in metadata and reports.

Before writing, check that no valid Float32 output elevation collides with the reserved `-255` value. If a collision occurs, stop with a diagnostic and revise the storage representation; never change a valid height to make room for the sentinel.

Readers must use declared nodata or an explicit validity mask. An `isfinite()` check alone is insufficient because `-255` is finite. NaN or masked arrays may still be used internally.

## River treatment and final COGs

Keep the coastal-only sheets as the working baseline. Use **GLO-30 WBM class 3**
for the river footprint, then derive water elevations from valid EGM2008-corrected
CartoDEM measurements around the river. [GLO30_INTEGRATION.md](river_flattening/GLO30_INTEGRATION.md)
records this source selection and its coverage checks. The shared-junction fit
developed during the WRIS pilots remains the fitting method.

Read categorical masks at native resolution, verifying their grid against each
CartoDEM sheet. Only class 3 forms river geometry; class 1 clips ocean and class 2
retains lake elevations without contributing bank observations. Class 0 supplies
eligible bank samples. Do not widen the mask, bridge land/islands, or treat
missing classification as land. WRIS is retained only as an explicit historical
comparison option, not mixed with GLO-30 at national borders.

The river stage should:

1. Identify connected reaches and obtain robust bank-height estimates along them.
2. Fit a non-increasing downstream water profile, with consistent elevations at junctions.
3. Make cross-sections level across the river, following that longitudinal profile.
4. Process reaches across sheet boundaries before writing results back into sheets.
5. Preserve unresolved nodata and report reaches with insufficient elevation or flow-direction information.
6. Fit connected sections with shared boundary elevations. Estimate those levels
   jointly from bank observations; do not pin every boundary pixel to its noisy
   original value. Keep an unresolved connected group unchanged if the model
   cannot treat its adjoining sections coherently. Record any direction inferred
   from the joint fit as a model hypothesis and audit the written raster joins.

Retain elevations outside GLO-30 river class 3. Export the river-treated result as a separate COG product so it can be compared with the coastal-only sheets. Do not substitute GLO-30 water heights. Lake class 2 is excluded from river treatment and bank samples.

Final COGs retain native Float32 elevations and explicitly declare **NoData=-255**. Both original source sentinels and ocean-masked pixels use that value; river processing preserves their validity. Use GDAL's COG writer, lossless DEFLATE with floating-point predictor, 512-pixel blocks and averaged numeric overviews. Verify every native elevation and validity sample against the river-treated working raster before accepting the COG. Only these final COGs are the publishable intermediate elevation product.

Method selection must follow published hydro-flattening work. [USGS surface-treatment guidance](https://www.usgs.gov/ngp-standards-and-specifications/lidar-base-specification-digital-elevation-model-surface) requires level cross-sections, downhill profiles and water edges no higher than nearby terrain. The [GRASS river workflow](https://grass-tutorials.osgeo.org/content/tutorials/hydro_flattening/hydro_flatenning.html#adding-break-lines) uses centreline transects and breaklines to subdivide the river before estimating elevations from its edges. Its lidar examples and percentile choices require validation on the coarser, noisier CartoDEM surface; they are not adopted as universal parameters. A river mask alone does not establish flow direction or a reliable longitudinal profile.

## Elevation resampling and tile coverage

Build a VRT mosaic of the selected intermediate sheets before tiling. This allows interpolation across sheet boundaries without storing a second complete elevation mosaic.

Use **512 × 512 XYZ tiles at zooms 0–12**, in Web Mercator, with a shared pixel grid at each zoom. This matches [Mapterhorn's global pyramid](https://mapterhorn.com/data-access/); its regional zooms 13–17 target higher-resolution sources. Use averaging of valid heights at zooms 0–11 and bilinear sampling at zoom 12. Generate every zoom from numeric elevations rather than resizing packed RGB tiles.

Carry a separate native 0/255 validity raster, with both values treated as data and no overviews. Derive it from the COG at native resolution. An implicit COG nodata mask can inherit averaged-overview behavior and hide small holes, so the tiler must not use that resampled implicit mask. Use nearest-neighbour sampling for the explicit validity raster and require both a valid mask sample and a finite resampled height for an output pixel to be valid. This is a pixel-centre coverage rule; features smaller than an output pixel may disappear at lower zooms. Check coastal strips, small islands and nodata boundaries in the pilot.

Nodata does not contribute `-255` or a zero-height placeholder to elevation averaging. Insert the zero-height placeholder only when encoding the final tile. Retain fully nodata tiles inside the selected tile extent as opaque zero-height tiles.

## Quantization and TerrainRGB encoding

Retain the z5–12 schedule in [the previous build commands](../steps.sh), and cap the new lower zooms at z5 precision:

| Zoom | Low bits cleared | Height increment in metres |
| --- | --- | --- |
| 0–4 | 11 | 204.8 |
| 5 | 11 | 204.8 |
| 6 | 10 | 102.4 |
| 7 | 9 | 51.2 |
| 8 | 8 | 25.6 |
| 9 | 7 | 12.8 |
| 10 | 6 | 6.4 |
| 11 | 5 | 3.2 |
| 12 | 4 | 1.6 |

For valid pixels, quantize the packed height code after elevation resampling, using `round_bits = 16 - max(zoom, 5)`. With nearest rounding, the added encoding error is at most half the listed increment, allowing for floating-point tolerance. The cap avoids the 6,553.6 m increment that extrapolating the old formula to z0 would produce. These increments describe storage precision, not source DEM accuracy.

Decode Mapbox TerrainRGB as:

```text
H = -10000 + 0.1 * (65536 * R + 256 * G + B)
```

| Pixel type | RGB | Alpha |
| --- | --- | --- |
| Valid terrain | Quantized height code | 255 |
| Nodata, including masked ocean | Exactly `(1, 134, 160)`, decoding to 0 m | 255 |

Apply the fixed nodata RGBA value **after quantization**, or exclude nodata from quantization entirely. It must remain `(1, 134, 160, 255)` at every zoom, even though those RGB bytes do not follow the cleared-bit pattern. Validate output against the explicit input validity mask; alpha no longer identifies valid pixels. A genuine zero-height terrain sample follows the normal quantization schedule.

Write **lossless WebP with every pixel opaque**. Check decoded RGBA equality after compression and after packaging. This supports normal TerrainRGB image loading. The earlier alpha-254 convention from [Iconem Terrain Viewer](https://terrain-viewer.iconem.com/docs/dev/demdiff-protocol/#holes-and-why-they-are-alpha-254) has been dropped; normal nodata remains in the COGs for analytical use.

**Browser validation, 25 October 2026:** default WebP decoding in Chrome 155
changes the red byte of `(1,134,160,254)` to zero, producing -6,553.6 m in
MapLibre 5.24.0. Do not assume ignoring alpha preserves RGB. The preview's
custom protocol uses `createImageBitmap` with `premultiplyAlpha: 'none'` and
`colorSpaceConversion: 'none'` and returns that bitmap directly to MapLibre.
That historical workaround preserves the original nodata marker and displays its placeholder at 0 m.
See the README's preview section and `preview/validation/` for the reproduction.
The new opaque display tiles do not require this custom decoder.

Package the tiles into PMTiles, retaining the image bytes. Metadata must identify Mapbox encoding, 512-pixel tile size, zooms 0–12, EGM2008 heights, attribution and the opaque zero-height nodata substitution. Associate each archive with its build manifest and quantization settings. For a completed z5–12 build with matching encoding and software, `extend_lowzooms.py` adds a separate z0–4 archive, verifies source and existing archive checksums, and atomically updates the mosaic manifest while retaining the original files and generation recipes.

## Validation and batch execution

Start with D43V around **12.2033 N, 75.1332 E**, plus neighbouring sheets to check seams. Include source samples containing each nodata sentinel and the previously inspected E45C and I41G locations. Evaluate the river variant separately after the coastal output is established.

Verify that:

- Both source sentinels are excluded before geoid correction, and valid zero and negative elevations remain eligible.
- Intermediate sheets retain native geometry, declare `-255` nodata, and contain no valid sentinel collisions.
- Ocean masking removes the intended pixels while preserving valid elevations elsewhere; land gaps remain identifiable.
- Reprojection excludes nodata from height interpolation and respects the selected coverage rule at coastlines and sheet edges.
- Every invalid output pixel is exactly `(1, 134, 160, 254)` at every zoom; every valid pixel has alpha 255.
- Valid decoded heights stay within the selected quantization error bound, and WebP and PMTiles preserve the encoded RGBA values.
- Supported river profiles are level across channels, non-increasing downstream and continuous at junctions and sheet boundaries; unsupported components retain their baseline and appear in reports.

Record timings, memory use and output sizes from the pilot to set worker counts for the full run. Use temporary output names followed by atomic completion, and checkpoint individual sheets and tile batches. Resume only outputs whose input hashes and processing settings match the manifest. Preserve intermediate sheets until the packaged output has passed validation.

Write run notes containing the exact commands, dependency versions, input inventory and checksums, processing settings, any deviations from this plan, validation results and output checksums. Keep connection details and machine-specific paths out of scripts and notes. A fresh working directory must be sufficient to reproduce the download and build; reusing verified existing downloads is an optimization.

## Implementation changes

| Component | Required work |
| --- | --- |
| [cartodem_to_egm2008.py](cartodem_to_egm2008.py) | Add mask-aware sheet preparation, ocean masking, uniform `-255` output and sentinel collision checks. Retain the existing geoid sampling and source validity handling. |
| [production_tiles.py](production_tiles.py) | Resample final numeric elevations and explicit validity separately, retain quantization, and support resumable tile batches. |
| [terrain_rgb.py](terrain_rgb.py) | Reuse its existing fixed nodata RGBA, quantization and lossless WebP behavior. Supply an explicit mask or masked elevations. |
| Batch driver and packaging | Add input manifests, bounded parallelism, completion records, reporting and PMTiles assembly. |
| [river_flattening](river_flattening/README.md) | Assemble geometry across sheets, fit shared CartoDEM-derived profiles, record unsupported components and apply supported surfaces. |
| [cog.py](cog.py) | Export final elevation COGs and verify native samples, nodata, overviews and COG layout. |
| [test_pipeline.py](test_pipeline.py) | Extend the existing checks for finite nodata, mask alignment, coastal masking, coverage through reprojection and packaged tile preservation. |

The README separates the final build commands from the historical standalone NaN-based converter examples. Keep the river method record and run notes current as validation proceeds.
