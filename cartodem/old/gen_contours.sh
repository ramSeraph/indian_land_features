#!/bin/bash

set -x

file=$1

tif_name=$(basename $1)
tname=${tif_name/%.tif}

echo $tname
[[ -e data/contours/${tname}.mbtiles && ! -e data/contours/${tname}.mbtiles-journal ]] && exit 0

#gdal_contour -snodata -255 -a ele_m -i 10 -f GeoJSONSeq $file data/contours/${tname}.geojsonl
gdal_contour -snodata -255 -a ele_m -i 10 -f FlatGeoBuf $file data/contours/${tname}.fgb

tippecanoe  -Z6  -z9 -P -y ele_m -l contour -C "jq 'select(.properties.ele_m % 500 == 0)'" data/contours/${tname}.fgb -o data/contours/${tname}-z6-9.mbtiles
tippecanoe -Z10 -z10 -P -y ele_m -l contour -C "jq 'select(.properties.ele_m % 200 == 0)'" data/contours/${tname}.fgb -o data/contours/${tname}-z10.mbtiles
tippecanoe -Z11 -z11 -P -y ele_m -l contour -C "jq 'select(.properties.ele_m % 100 == 0)'" data/contours/${tname}.fgb -o data/contours/${tname}-z11.mbtiles
tippecanoe -Z12 -z12 -P -y ele_m -l contour -C "jq 'select(.properties.ele_m %  50 == 0)'" data/contours/${tname}.fgb -o data/contours/${tname}-z12.mbtiles
tippecanoe -Z13 -z13 -P -y ele_m -l contour -C "jq 'select(.properties.ele_m %  20 == 0)'" data/contours/${tname}.fgb -o data/contours/${tname}-z13.mbtiles
tippecanoe -Z14 -z14 -P -y ele_m -l contour data/contours/${tname}.fgb -o data/contours/${tname}-z14.mbtiles
tile-join -o data/contours/${tname}.mbtiles data/contours/${tname}-z*.mbtiles

rm data/contours/${tname}-z*.mbtiles
rm data/contours/${tname}.fgb
