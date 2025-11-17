#!/bin/bash

uvx --from wmsdump wms-extractor extract -u https://uat.soiindia.in/geoserver/soi/wms soi:CONTOUR_Crv -s wms -b 1000 -m extent --bounds "67.38,7.97,97.37,35.37" --max-box-dims "0.1,0.1"

uvx --from wmsdump geojsonl-dedupe soi_CONTOUR_Crv.geojsonl

mv deduped_soi_CONTOUR_Crv.geojsonl SOI_Contours_raw.geojsonl

./tile.sh

pmtiles convert SOI_Contours_raw.mbtiles SOI_Contours_raw.pmtiles


