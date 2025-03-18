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

uv venv
uv pip install "git+https://github.com/acalcutt/rio-rgbify"
uv run rio rgbify -b -10000 -i 0.1 --min-z  5 --max-z  5 --round-digits 11 -j 8 --format webp data/dem_epsg3857.vrt data/cartodem_v3r1_z5.mbtiles
uv run rio rgbify -b -10000 -i 0.1 --min-z  6 --max-z  6 --round-digits 10 -j 8 --format webp data/dem_epsg3857.vrt data/cartodem_v3r1_z6.mbtiles
uv run rio rgbify -b -10000 -i 0.1 --min-z  7 --max-z  7 --round-digits  9 -j 8 --format webp data/dem_epsg3857.vrt data/cartodem_v3r1_z7.mbtiles
uv run rio rgbify -b -10000 -i 0.1 --min-z  8 --max-z  8 --round-digits  8 -j 8 --format webp data/dem_epsg3857.vrt data/cartodem_v3r1_z8.mbtiles
uv run rio rgbify -b -10000 -i 0.1 --min-z  9 --max-z  9 --round-digits  7 -j 8 --format webp data/dem_epsg3857.vrt data/cartodem_v3r1_z9.mbtiles
uv run rio rgbify -b -10000 -i 0.1 --min-z 10 --max-z 10 --round-digits  6 -j 8 --format webp data/dem_epsg3857.vrt data/cartodem_v3r1_z10.mbtiles
uv run rio rgbify -b -10000 -i 0.1 --min-z 11 --max-z 11 --round-digits  5 -j 8 --format webp data/dem_epsg3857.vrt data/cartodem_v3r1_z11.mbtiles
uv run rio rgbify -b -10000 -i 0.1 --min-z 12 --max-z 12 --round-digits  4 -j 8 --format webp data/dem_epsg3857.vrt data/cartodem_v3r1_z12.mbtiles


uvx --from "git+https://github.com/mapbox/mbutil" mb-util --image_format webp data/cartodem_v3r1_z5.mbtiles data/tiles5
uvx --from "git+https://github.com/mapbox/mbutil" mb-util --image_format webp data/cartodem_v3r1_z6.mbtiles data/tiles6
uvx --from "git+https://github.com/mapbox/mbutil" mb-util --image_format webp data/cartodem_v3r1_z7.mbtiles data/tiles7
uvx --from "git+https://github.com/mapbox/mbutil" mb-util --image_format webp data/cartodem_v3r1_z8.mbtiles data/tiles8
uvx --from "git+https://github.com/mapbox/mbutil" mb-util --image_format webp data/cartodem_v3r1_z9.mbtiles data/tiles9
uvx --from "git+https://github.com/mapbox/mbutil" mb-util --image_format webp data/cartodem_v3r1_z10.mbtiles data/tiles10
uvx --from "git+https://github.com/mapbox/mbutil" mb-util --image_format webp data/cartodem_v3r1_z11.mbtiles data/tiles11
uvx --from "git+https://github.com/mapbox/mbutil" mb-util --image_format webp data/cartodem_v3r1_z12.mbtiles data/tiles12

mkdir data/tiles
cp data/tiles5/metadata.json data/tiles/

mv data/tiles5/5 data/tiles
mv data/tiles6/6 data/tiles
mv data/tiles7/7 data/tiles
mv data/tiles8/8 data/tiles
mv data/tiles9/9 data/tiles
mv data/tiles10/10 data/tiles
mv data/tiles11/11 data/tiles
mv data/tiles12/12 data/tiles

rm -rf data/tiles5
rm -rf data/tiles6
rm -rf data/tiles7
rm -rf data/tiles8
rm -rf data/tiles9
rm -rf data/tiles10
rm -rf data/tiles11
rm -rf data/tiles12

uvx --from "git+https://github.com/mapbox/mbutil" mb-util --image_format webp data/tiles Bhuvan_CartoDEM_v3r1_TerrainRGB.mbtiles
