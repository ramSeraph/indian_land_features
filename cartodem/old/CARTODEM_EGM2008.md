# Converting original CartoDEM elevations to EGM2008 for tile generation

Original CartoDEM elevations can be converted to EGM2008 heights by subtracting a sampled EGM2008 geoid undulation from each valid source pixel: **H = h − N**. Preserve the original raster grid, store the result as Float32 metres, and carry validity explicitly through later reprojection and tile generation.

This procedure uses original CartoDEM TIFFs and one geoid grid. Copernicus GLO-30 is an optional validation reference. The accompanying converter preserves CartoDEM coverage and uses only its own validity information.

## Inputs and datum assumption

Use the original rasters from the [CartoDEM 30 m V3R1 original release](https://github.com/ramSeraph/indian_land_features/releases/tag/cartodem-30m-v3r1-orig). Its `listing_files.csv` lists the assets; `grid.geojson` provides the sheet index. Locally, the source directory is `original_cartodem/`, containing 543 TIFFs.

The four inspected sheets are single-band Int16 rasters in EPSG:4326, with 3600 × 3600 pixels at one arcsecond spacing. Their TIFFs do not declare a vertical CRS or a nodata value. Integer samples have 1 m numeric increments.

**Source-height assumption:** interpret the original values as WGS84 ellipsoidal heights in metres. Table 23 of the [Bhuvan Geospatial Content Standards](https://bhuvan-app1.nrsc.gov.in/2dresources/documents/4_Bhuvan_Data_Content_And_Map_Standards.pdf), printed page 43 / PDF page 47, lists WGS84 for the CartoDEM vertical datum. The independent comparisons below support this interpretation for these sheets; it remains an explicit assumption for a release without per-file vertical metadata.

The target is EGM2008 height, identified by vertical CRS EPSG:3855. EPSG:4326 describes the horizontal coordinates; assigning it does not perform a height conversion. Apply this procedure once to original elevations. An input already expressed in EGM2008 heights must not be corrected again.

## The geoid grid and correction sign

Use [us_nga_egm08_25.tif](https://cdn.proj.org/us_nga_egm08_25.tif), the 2.5 arcminute EGM2008 grid distributed by PROJ. It is a geoid undulation grid derived from NGA data, with provenance described in the [PROJ grid README](https://cdn.proj.org/us_nga_README.txt).

The exact file used here is 80,585,622 bytes, with SHA-256:

```text
4191d471eefebf24091b56dbc604353cb3b8cf8cc70e448bb9ae56a272bef17a
```

Definitions, all in metres:

- `h`: original CartoDEM height, assumed WGS84 ellipsoidal.
- `N`: EGM2008 geoid undulation sampled at that pixel.
- `H`: resulting EGM2008 height.

```text
h = H + N
H = h - N
```

In the four checked sheets, N is negative, so subtracting it increases the height. For example, if h = −66 m and N = −67 m, H = +1 m. This is an illustration of the sign convention, not a constant offset to apply everywhere.

| Sheet | Sampled N over valid CartoDEM pixels | Height added by subtracting N |
| --- | --- | --- |
| Delhi h43x | −56.68 to −51.15 m | +51.15 to +56.68 m |
| Leh i43l | −25.37 to −22.68 m | +22.68 to +25.37 m |
| Mumbai e43a | −68.72 to −65.09 m | +65.09 to +68.72 m |
| Chennai d44o | −92.27 to −90.14 m | +90.14 to +92.27 m |

Sample N spatially for every pixel; a sheet-wide mean or a single India-wide offset will introduce errors. The geoid grid is 150 times coarser in angular spacing than these one arcsecond DEMs. Its interpolation adds a smooth height correction; it does not add terrain detail or improve the source's vertical accuracy.

## Exact raster procedure

The geoid subtraction in the original comparison code, `compare_glo30.py`, used the following procedure. The standalone `cartodem_to_egm2008.py` carries out the same subtraction and sampling with CartoDEM's own validity mask.

1. Read the original TIFF's CRS, affine transform, width and height. Preserve them exactly.
2. Interpret the geoid grid as a scalar raster of N values on a two-dimensional longitude–latitude grid.
3. Bilinearly sample N onto the original CartoDEM grid, at the source pixel centres.
4. Read original heights as Float32, identify valid pixels, and compute `H = h - N` there.
5. Write Float32 GeoTIFFs in metres, with NaN nodata, 256 × 256 storage blocks, DEFLATE compression and floating-point predictor 3.

The key sampling configuration is:

```python
with WarpedVRT(
    geo,
    src_crs="EPSG:4326",
    crs=src.crs,
    transform=src.transform,
    width=src.width,
    height=src.height,
    dtype="float32",
    resampling=Resampling.bilinear,
    nodata=np.nan,
    warp_mem_limit=16,
    tolerance=1e-6,
) as geoid:
    n = geoid.read(1, window=window)
    h = src.read(1, window=window).astype(np.float32)
    corrected = h - n
```

The actual conversion applies the validity mask before writing corrected values. See the complete script below.

**Why override the geoid source CRS?** The grid carries EPSG:4979 metadata, a three-dimensional WGS84 CRS. Here its pixel values must be sampled as scalar undulations, followed by one explicit subtraction. The two-dimensional `src_crs="EPSG:4326"` override prevents this sampling step from invoking a separate vertical transformation. It changes the reader's interpretation for this operation; it does not edit the geoid file. [Rasterio WarpedVRT parameters](https://rasterio.readthedocs.io/en/stable/api/rasterio.vrt.html), [GDAL vertical transformation behavior](https://gdal.org/en/stable/programs/gdalwarp.html#vertical-transformation).

The CartoDEM height raster itself is not resampled during correction. Preserve its numeric transform, including its existing half-pixel offsets. Do not recreate the transform from integer sheet bounds or shift the raster because another product uses different `AREA_OR_POINT` metadata.

Processing uses 512 × 512 windows, a 24 MiB GDAL cache, a 16 MB warp buffer and one GDAL thread. The output retains horizontal EPSG:4326 and records the vertical target and source assumption in metadata tags. Those tags describe the height values; they do not constitute a formal compound CRS.

## Nodata and geographic coverage

For this original release, **−255 is treated as inferred nodata**, despite the TIFFs lacking a declared nodata value. Large coastal fill regions in the inspected originals support that interpretation. Confirm the convention if using a different CartoDEM release.

The conversion mask is:

```python
source_valid = (
    (src.read_masks(1, window=window) != 0)
    & np.isfinite(h)
    & (h != -255)
)
valid = source_valid & np.isfinite(n)
```

Identify −255 in the original values before correction. Subtracting N first would turn fill into ordinary-looking negative elevations and lose the exact sentinel.

Zero and other negative source elevations remain eligible. Negative ellipsoidal heights are expected in some low coastal areas when N is negative. Corrected values below zero also remain unchanged; datum conversion does not decide whether an individual terrain measurement is reliable.

The converter writes invalid pixels as NaN and declares NaN nodata. It does not fill holes, flatten water or clip coverage. If an India-only product is desired, apply the [DataMeet India boundary](https://github.com/datameet/maps/blob/master/Country/india-soi.geojson) as a separate coverage policy and record the boundary version and pixel-inclusion rule. That policy is independent of the height datum.

For tools that cannot preserve NaN reliably, export an explicitly declared finite nodata sentinel outside the valid height range, or carry a separate validity mask. Ensure every subsequent operation honors that representation.

## Run the conversion elsewhere

Save the complete script below as `cartodem_to_egm2008.py`. The only Python packages it needs are NumPy and Rasterio. The verification run used NumPy 2.5.4, Rasterio 1.5.2, GDAL 3.12.2 and PROJ 9.8.1.

For example, on a machine with Python and curl:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install numpy==2.5.4 rasterio==1.5.2

mkdir -p original_cartodem geoid corrected_cartodem

curl -fL --retry 3 \
  https://cdn.proj.org/us_nga_egm08_25.tif \
  -o geoid/us_nga_egm08_25.tif

curl -fL --retry 3 \
  https://github.com/ramSeraph/indian_land_features/releases/download/cartodem-30m-v3r1-orig/cdne43a.tif \
  -o original_cartodem/cdne43a.tif

.venv/bin/python cartodem_to_egm2008.py \
  original_cartodem/cdne43a.tif \
  geoid/us_nga_egm08_25.tif \
  corrected_cartodem/cdne43a_egm2008.tif
```

The script verifies the pinned geoid checksum, refuses to overwrite an output, checks the original raster format, and records source and grid hashes with runtime versions. It fails if geoid samples are unavailable at otherwise-valid input pixels.

For multiple downloaded original sheets, process them sequentially to keep memory bounded:

```bash
for src in original_cartodem/*.tif; do
  name="${src##*/}"
  stem="${name%.tif}"
  .venv/bin/python cartodem_to_egm2008.py \
    "$src" geoid/us_nga_egm08_25.tif \
    "corrected_cartodem/${stem}_egm2008.tif" || break
done
```

Keep the original directory separate from corrected outputs. Runtime GDAL and PROJ versions can differ across platform-specific Rasterio installations; the output tags and JSON report record what actually ran.

## Independent validation

The standalone implementation was checked on all pixels of four complete source sheets on 21 October 2026. The comparison used only valid original CartoDEM pixels, finite geoid samples and valid GLO-30 pixels. There was no country-boundary filter.

GLO-30 is already referenced to EGM2008, so it received no height correction. For these four downloaded pairs, the numeric dimensions, CRS and affine transforms matched exactly. That must be checked again for other files. [Copernicus DEM specifications](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM).

| Sheet | Compared pixels | Raw RMSE vs GLO-30 | Corrected RMSE vs GLO-30 | Corrected mean difference |
| --- | ---: | ---: | ---: | ---: |
| Delhi h43x | 12,960,000 | 52.66 m | 2.87 m | +1.21 m |
| Leh i43l | 12,960,000 | 86.88 m | 80.92 m | −8.37 m |
| Mumbai e43a | 3,197,060 | 67.71 m | 5.55 m | −0.70 m |
| Chennai d44o | 3,571,066 | 89.81 m | 3.00 m | +1.96 m |

Differences are CartoDEM minus GLO-30. These results support the correction's sign and magnitude, particularly on the lowland sheets. Leh retains large residual errors, showing that a datum correction cannot resolve all differences in mountainous terrain. Both products are surface models; GLO-30 is a comparison reference, not ground truth.

Verification also established:

- 32,688,126 source-valid pixels remained finite after correction; all 19,151,874 inferred source-fill pixels remained invalid.
- All 21,718 raw zero-valued pixels in Mumbai and Chennai remained eligible for correction.
- Independent interpolation directly from the geoid's four neighboring cells agreed with 80 output samples to within 0.000244 m, consistent with Float32 rounding.
- Raster geometry and `AREA_OR_POINT` were preserved. No otherwise-valid source pixel lacked a geoid sample.
- Peak resident memory for the validation process was 142.7 MiB.

The machine-readable record is `analysis/cartodem_egm2008_standalone_validation.json`. It includes source hashes, sampled coordinates and heights, nodata counts, version numbers and comparison metrics. `validate_cartodem_egm2008.py` reproduces the checks using temporary output rasters.

## Prepare corrected elevations for tiles

Keep the corrected Float32 elevation rasters as the authoritative intermediate. A suitable processing order is:

```text
Original CartoDEM with explicit validity
    -> EGM2008 correction on the original grid
    -> mosaic of corrected elevations
    -> numeric elevation reprojection and pyramid
    -> encode each zoom into terrain RGB
    -> lossless image compression
```

Build the mosaic before resampling tiles so interpolation can use adjacent sheets. Sample a halo around tile boundaries when computing hillshade, slope or contours, then crop to the tile footprint. Preserve nodata throughout.

For an elevation pyramid, resample numeric elevations. Bilinear interpolation is a reasonable starting point for reprojection, and an average of valid elevations is a reasonable starting point for coarser levels. Keep a validity or coverage raster alongside the elevations. If conservative edge coverage is required, resample a separate 0/1 validity raster with zero treated as data, then require its support fraction to meet a documented threshold.

After the explicit H = h − N conversion, later reprojection must not apply the geoid shift again. For a GDAL horizontal warp, `-novshift` disables automatic vertical transformation. EPSG:3857 specifies the horizontal web-map projection; it does not replace the EGM2008 height reference. [GDAL gdalwarp documentation](https://gdal.org/en/stable/programs/gdalwarp.html#cmdoption-gdalwarp-novshift).

For Mapbox-compatible Terrain-RGB, encode each valid height at the final zoom resolution:

```python
# H is the numeric EGM2008 elevation raster; valid is its boolean mask.
q = np.rint((H[valid].astype(np.float64) + 10000.0) * 10.0)
if np.any((q < 0) | (q > 16777215)):
    raise ValueError("Elevation is outside the Terrain-RGB encoding range")
q = q.astype(np.uint32)
R = ((q >> 16) & 255).astype(np.uint8)
G = ((q >> 8) & 255).astype(np.uint8)
B = (q & 255).astype(np.uint8)
```

Decoding is `H = -10000 + 0.1 * (65536*R + 256*G + B)`, following the [Terrain-RGB encoding specification](https://docs.mapbox.com/data/tilesets/reference/mapbox-terrain-rgb-v1/). Rounding adds at most about 0.05 m of encoding error for valid heights. This is storage precision, not DEM accuracy.

Do not generate lower zooms by resizing encoded RGB images. Each channel is a digit of a packed height value, and image resampling rounds those digits independently. Generate the numeric elevation for each zoom first, then encode it. Avoid image color correction, JPEG, lossy WebP and palette quantization. For lossless WebP, verify decoded RGB equality on valid pixels.

Terrain-RGB's height formula has no inherent nodata flag. Select a validity representation that the intended renderer actually honors, such as a supported alpha mask or a separate mask layer. A zero RGB triplet decodes to −10,000 m; it is not automatically missing data. A height of zero metres is also valid data. Test partially valid tiles and fully empty tiles explicitly.

For cartographic output, keep height-datum conversion, void filling, water treatment, smoothing and contour generalization as separate recorded steps. The datum correction addresses the reference surface; it does not remove surface-model artifacts or make an unedited DEM suitable for every map style.

## Portable conversion script

Save this block as `cartodem_to_egm2008.py`.

```python
#!/usr/bin/env python3
"""Convert original CartoDEM Int16 sheets from assumed ellipsoidal heights to EGM2008.

Usage:
  python cartodem_to_egm2008.py INPUT.tif GEOID.tif OUTPUT.tif
Uses only CartoDEM and the pinned EGM2008 grid. No comparison or boundary mask.
"""
import os
os.environ.setdefault("GDAL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse
import hashlib
import json
import tempfile
from pathlib import Path
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window

GEOID_SHA256 = "4191d471eefebf24091b56dbc604353cb3b8cf8cc70e448bb9ae56a272bef17a"

def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def windows(width, height):
    for y in range(0, height, 512):
        for x in range(0, width, 512):
            yield Window(x, y, min(512, width-x), min(512, height-y))

def convert(source_path, geoid_path, output_path, source_nodata=-255):
    source_path, geoid_path, output_path = map(Path, (source_path, geoid_path, output_path))
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite {output_path}")
    if output_path.resolve() in (source_path.resolve(), geoid_path.resolve()):
        raise ValueError("Output must differ from both inputs")
    geoid_hash = sha256(geoid_path)
    if geoid_hash != GEOID_SHA256:
        raise ValueError("Geoid differs from the pinned us_nga_egm08_25.tif; verify its model, units and sign")
    source_hash = sha256(source_path)
    temporary = None
    valid_count = nodata_count = geoid_missing = 0
    nmin, nmax = np.inf, -np.inf
    hmin, hmax = np.inf, -np.inf
    try:
        with rasterio.Env(GDAL_CACHEMAX=24*1024**2, GDAL_NUM_THREADS="1"), \
             rasterio.open(source_path) as src, rasterio.open(geoid_path) as geo:
            if src.count != 1 or src.dtypes != ("int16",) or src.crs != rasterio.crs.CRS.from_epsg(4326):
                raise ValueError("Expected an original, single-band Int16 CartoDEM sheet in EPSG:4326")
            if src.scales != (1.0,) or src.offsets != (0.0,):
                raise ValueError("Unexpected scale/offset; inspect the input before conversion")
            if src.transform.b != 0 or src.transform.d != 0:
                raise ValueError("Expected a north-up original sheet")
            if "HEIGHT_FORMULA" in src.tags() or "VERTICAL_DATUM" in src.tags():
                raise ValueError("Input already has height-datum processing tags; inspect before converting")
            profile = dict(driver="GTiff", width=src.width, height=src.height, count=1,
                           dtype="float32", crs=src.crs, transform=src.transform,
                           nodata=np.nan, tiled=True, blockxsize=256, blockysize=256,
                           compress="DEFLATE", predictor=3, BIGTIFF="IF_SAFER")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(prefix=output_path.stem+".", suffix=".part.tif",
                                        dir=output_path.parent)
            os.close(fd)
            temporary = Path(name)
            # Sample scalar N only. The grid's EPSG:4979 metadata must not
            # trigger an automatic vertical transform during this resampling.
            with WarpedVRT(geo, src_crs="EPSG:4326", crs=src.crs,
                           transform=src.transform, width=src.width, height=src.height,
                           dtype="float32", resampling=Resampling.bilinear,
                           nodata=np.nan, warp_mem_limit=16, tolerance=1e-6) as geoid, \
                 rasterio.open(temporary, "w", **profile) as dst:
                dst.set_band_description(1, "EGM2008 height under WGS84 ellipsoidal source assumption")
                dst.set_band_unit(1, "m")
                if "AREA_OR_POINT" in src.tags():
                    dst.update_tags(AREA_OR_POINT=src.tags()["AREA_OR_POINT"])
                dst.update_tags(
                    HEIGHT_FORMULA="H = h - N",
                    SOURCE_VERTICAL_DATUM="Assumed WGS84 ellipsoidal heights in metres",
                    VERTICAL_DATUM="EGM2008 under stated source assumption",
                    TARGET_VERTICAL_CRS="EPSG:3855",
                    GEOID_MODEL="EGM2008",
                    GEOID_GRID="us_nga_egm08_25.tif",
                    GEOID_SHA256=geoid_hash,
                    GEOID_RESAMPLING="bilinear at original CartoDEM pixel centres",
                    SOURCE_FILE=source_path.name,
                    SOURCE_SHA256=source_hash,
                    SOURCE_NODATA_POLICY=f"existing mask or nonfinite or original value == {source_nodata}",
                    RASTERIO_VERSION=rasterio.__version__,
                    GDAL_VERSION=rasterio.__gdal_version__,
                    PROJ_VERSION=rasterio.__proj_version__,
                )
                for window in windows(src.width, src.height):
                    h = src.read(1, window=window).astype(np.float32)
                    n = geoid.read(1, window=window)
                    source_valid = (src.read_masks(1, window=window) != 0) & np.isfinite(h) & (h != source_nodata)
                    valid = source_valid & np.isfinite(n)
                    result = np.full(h.shape, np.nan, dtype=np.float32)
                    np.subtract(h, n, out=result, where=valid)
                    dst.write(result, 1, window=window)
                    valid_count += int(valid.sum())
                    nodata_count += int((~source_valid).sum())
                    geoid_missing += int((source_valid & ~np.isfinite(n)).sum())
                    if valid.any():
                        nmin = min(nmin, float(n[valid].min()))
                        nmax = max(nmax, float(n[valid].max()))
                        hmin = min(hmin, float(result[valid].min()))
                        hmax = max(hmax, float(result[valid].max()))
            if geoid_missing:
                raise ValueError(f"Geoid unavailable at {geoid_missing} otherwise-valid input pixels")
            with rasterio.open(temporary) as check:
                assert check.transform == src.transform and check.crs == src.crs and check.shape == src.shape
                assert check.dtypes == ("float32",) and np.isnan(check.nodata)
                written_valid = sum(int(np.isfinite(check.read(1, window=w)).sum())
                                    for w in windows(check.width, check.height))
                assert written_valid == valid_count
            if output_path.exists():
                raise FileExistsError(f"Output appeared during conversion: {output_path}")
            temporary.replace(output_path)
            temporary = None
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return dict(source=str(source_path), output=str(output_path),
                source_sha256=source_hash, geoid_sha256=geoid_hash,
                valid_pixels=valid_count, source_nodata_pixels=nodata_count,
                geoid_missing_pixels=geoid_missing,
                geoid_undulation_range_m=[nmin,nmax] if valid_count else None,
                corrected_height_range_m=[hmin,hmax] if valid_count else None,
                formula="H = h - N",
                vertical_assumption="Original CartoDEM is WGS84 ellipsoidal; target EGM2008",
                versions=dict(numpy=np.__version__, rasterio=rasterio.__version__,
                              gdal=rasterio.__gdal_version__, proj=rasterio.__proj_version__))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("geoid", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--source-nodata", type=float, default=-255,
                        help="Inferred fill in this original CartoDEM release (default: -255)")
    args = parser.parse_args()
    print(json.dumps(convert(args.source, args.geoid, args.output, args.source_nodata),
                     indent=2, allow_nan=False))

if __name__ == "__main__":
    main()
```
