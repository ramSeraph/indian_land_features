# Derived data from [CartoDEM V3R1](https://github.com/ramSeraph/indian_land_features/releases/tag/cartodem-30m-v3r1-orig)

The COG rasters and TerrainRGB tiles in this release are derived from
**CartoDEM V3R1**, with heights converted to **EGM2008**, followed by ocean masking
using the **Copernicus GLO-30 ocean mask** and **experimental river flattening
using the GLO-30 river mask**. GLO-30 supplies the masks; river elevations are
estimated from the corrected CartoDEM data. CartoDEM's existing lake treatment
is retained.

River flattening is experimental. Unsupported rivers retain their corrected
CartoDEM elevations, and the [river method and measured limitations](egm2008/river_flattening/README.md)
document the remaining gradient issues.

## Cloud Optimized GeoTIFFs

The release contains 543 native-resolution Float32 COGs. They preserve missing
data as `NoData=-255`, including original `-1024` and `-255` pixels and pixels
removed by ocean masking.

Download filenames, sizes and links are in
[listing_files.csv](https://github.com/ramSeraph/indian_land_features/releases/download/cartodem-30m-v3r1/listing_files.csv).

## TerrainRGB tiles

The four PMTiles archives contain 512 × 512 lossless WebP tiles at **zooms 0–12**,
using Mapbox TerrainRGB encoding and zoom-dependent height quantization.
Missing data is displayed as opaque **0 m**, without an alpha-254 marker.
Use the COGs when the distinction between zero elevation and nodata matters.

- Tile URL: `https://indianopenmaps.com/dem/terrain-rgb/cartodem-v3r1/bhuvan/{z}/{x}/{y}.webp`
- [View the tiles](https://indianopenmaps.com/dem/terrain-rgb/cartodem-v3r1/bhuvan/rasterview)
- [PMTiles mosaic](https://github.com/ramSeraph/indian_land_features/releases/download/cartodem-30m-v3r1/Bhuvan_CartoDEM_v3r1_TerrainRGB.mosaic.json)

To download the mosaic as a local [MBTiles](https://docs.mapbox.com/help/glossary/mbtiles/)
file, install [uv](https://docs.astral.sh/uv/getting-started/installation) and run:

```bash
uvx --from pmtiles-mosaic download-mosaic \
  -u https://github.com/ramSeraph/indian_land_features/releases/download/cartodem-30m-v3r1/Bhuvan_CartoDEM_v3r1_TerrainRGB.mosaic.json \
  -o terrain_rgb.mbtiles
```

## Reproduce the data

The [scripts and reproduction instructions](egm2008/README.md) cover downloading
the inputs, EGM2008 conversion, ocean masking, experimental river fitting, COG
export, TerrainRGB encoding and partitioning with `pmtiles-mosaic`.
