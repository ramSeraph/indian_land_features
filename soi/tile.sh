#!/bin/bash

set -ex

fname=SOI_Contours_raw

tippecanoe  -Z6  -z10 -P -y VALUE -l contour -C "jq 'select((.properties.VALUE | tonumber | floor ) % 400 == 0)'" ${fname}.geojsonl -o ${fname}-z6-10.mbtiles
tippecanoe -Z10 -z10 -P -y VALUE -l contour -C "jq 'select((.properties.VALUE | tonumber | floor ) % 200 == 0)'" ${fname}.geojsonl -o ${fname}-z10.mbtiles
tippecanoe -Z11 -z11 -P -y VALUE -l contour -C "jq 'select((.properties.VALUE | tonumber | floor ) %  80 == 0)'" ${fname}.geojsonl -o ${fname}-z11.mbtiles
tippecanoe -Z12 -z12 -P -y VALUE -l contour -C "jq 'select((.properties.VALUE | tonumber | floor ) %  40 == 0)'" ${fname}.geojsonl -o ${fname}-z12.mbtiles
tippecanoe -Z13 -z13 -P -y VALUE -l contour -C "jq 'select((.properties.VALUE | tonumber | floor ) %  20 == 0)'" ${fname}.geojsonl -o ${fname}-z13.mbtiles
tippecanoe -Z14 -z14 -P -y VALUE -l contour ${fname}.geojsonl -o ${fname}-z14.mbtiles
tile-join -o ${fname}.mbtiles ${fname}-z*.mbtiles
