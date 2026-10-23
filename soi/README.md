# SOI contour-derived DEM — WIP

Intermediate elevation data derived automatically from the released Survey of
India contours. These sheets will be replaced by final data. Waterbody and river
flattening have not been applied.

- [Production workflow and reproduction instructions](PRODUCTION.md)
- [Historical calibration notes](calibration/NOTES.md)
- [Official statement on DEM and contour availability](PRODUCTION.md#official-data-availability-statement)
- [Intermediate elevation sheets](https://github.com/ramSeraph/indian_land_features/releases/tag/soi-dem-30m)

The workflow reads local contour GeoParquet partitions through iomaps, includes
all accepted height points, interpolates at 30 m, reconciles available shared
sheet boundaries, and exports Cloud Optimized GeoTIFFs. Circle and triangle
height points are retained separately as GeoJSONL. Toposheets and Copernicus
GLO-30 are used for calibration and comparison only.

The Python interpolation port is based on Christoph Hormann's
[Imagico contour interpolation](https://www.imagico.de/pov/earth_interpolate.php).
Its numerical implementation and raster adapter use GPL-3.0-or-later; see
[the license](LICENSE.imagico.txt) and [provenance](imagico_provenance.json).
Original source is retained in `contour_interpolate-0.3/`. The source data's
licensing statement is copied in [DATA_LICENSE.source.md](DATA_LICENSE.source.md).
Credit Survey of India and the Datameet community for the source data.
