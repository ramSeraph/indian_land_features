#!/bin/bash

# get the grid file
uvx --from wmsdump wms-extractor extract cdem_v3_r1 data/cdem_v3_r1_grid.geojsonl -u https://bhuvan-vec2.nrsc.gov.in/bhuvan/vector/wms 

# create data/user.json file with the following format
#{
#  "name": "<username>",
#  "password": "<password>"
#}

# scrape the data to data/raw/v3_r1
mkdir =p data/raw/v3_r1
uv run scrape_bhuvan.py
# you will be prompted to solve the captcha.. not worth automating at this point

# extract the tif file to data/tifs/v3_r1
mkdir -p data/tifs/v3_r1
cd data/tifs/v3_r1
ls ../../raw/v3_r1/ | cut -d"." -f1 | xargs -I {} 7z e ../../raw/v3_r1/{}.zip cdn{}_v3r1/cdn{}.tif
cd -

# collect and zip the water bodies data
cd data/raw/v3_r1
mkdir -p data/water/v3_r1
ls | cut -d"." -f1 | xargs -I {} ogr2ogr -f GeoJSONSeq ../../water/v3_r1/{}.geojsonl /vsizip/{}.zip/cdn{}_v3r1/cdn{}.shp
cd -
uvrun join.py
rm -f data/water/v3_r1/*
cd data
7z a v3_r1_water_bodies.geojsonl.7z v3_r1_water_bodies.geojsonl
rm v3_r1_water_bodies.geojsonl.7z
cd -

# create contour vector tiles
find data/tifs/v3_r1 -type f | xargs -I {} ./gen_contours.sh {}
tile-join -n Bhuvan_CartoDEM_v3r1_Contours -l contour -pk -o Bhuvan_CartoDEM_v3r1_Contours.mbtiles data/contours/*.mbtiles

# create terrain rgb tiles
gdalbuildvrt -vrtnodata -255 data/dem_255.vrt data/raw/v3_r1/*.tif
gdalwarp -r cubicspline -s_srs EPSG:4326 -t_srs EPSG:3857 -dstnodata 0 -co COMPRESS=DEFLATE data/dem_255.vrt data/dem_epsg3857.vrt


mkdir data/tiles
cp metadata.json data/tiles/
round_bits=11
for zoom in 5 6 7 8 9 10 11 12
do
    uvx --from rasterio --with "git+https://github.com/acalcutt/rio-rgbify" rio rgbify -b -10000 -i 0.1 --min-z  ${zoom} --max-z  ${zoom}  --round-digits $round_bits -j 8 --format webp data/dem_epsg3857.vrt data/cartodem_v3r1_z${zoom}.mbtiles

    uvx --from "git+https://github.com/mapbox/mbutil" mb-util --image_format webp data/cartodem_v3r1_z5.mbtiles data/tiles${zoom}

    mv data/tiles${zoom}/${zoom} data/tiles
    rm -rf data/tiles${zoom}

    round_bits=$(( round_bits - 1))
done

uvx --from "git+https://github.com/mapbox/mbutil" mb-util --image_format webp data/tiles Bhuvan_CartoDEM_v3r1_TerrainRGB.mbtiles
