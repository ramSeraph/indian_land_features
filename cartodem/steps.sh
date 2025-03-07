#!/bin/bash

uvx --from wmsdump wms-extractor extract cdem_v3_r1 data/cdem_v3_r1_grid.geojsonl -u https://bhuvan-vec2.nrsc.gov.in/bhuvan/vector/wms 

# create data/user.json file with the following format
#{
#  "name": "<username>",
#  "password": "<password>"
#}

uv run scrape_bhuvan.py
# you will be prompted to solve the captcha.. not worth automating at this point
