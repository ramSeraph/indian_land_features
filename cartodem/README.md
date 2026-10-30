# CartoDEM processing

[egm2008/](egm2008/README.md) contains the current reproducible workflow:

1. Download and verify original CartoDEM V3R1 sheets, the EGM2008 geoid grid,
   and Copernicus GLO-30 water classification rasters.
2. Unify nodata, correct heights to EGM2008, and remove GLO-30 ocean pixels.
3. Apply experimental river flattening using the GLO-30 river mask and
   CartoDEM bank elevations, with connected fits across sheet boundaries.
4. Export native-resolution Float32 COGs with `NoData=-255`.
5. Generate opaque TerrainRGB WebP tiles, export standard MBTiles, and use
   `pmtiles-mosaic` to create geographic PMTiles partitions.

The [river method](egm2008/river_flattening/README.md) explains fitting and its
limitations. [old/](old/README.md) preserves superseded scripts and experiments.
[RELEASE_NOTES.md](RELEASE_NOTES.md) contains the derived-data release notes.
