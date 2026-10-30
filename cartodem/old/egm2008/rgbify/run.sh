#!/usr/bin/env bash
# Corrected COGs -> GDAL VRT -> rio-rgbify MBTiles -> pmtiles-mosaic partitions.
# Uses the unmodified rio-rgbify quantization behavior; see README.md.
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 COG_DIRECTORY OUTPUT_DIRECTORY [WORKERS=4]" >&2
    exit 2
fi
script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
cog_dir=$(cd -- "$1" && pwd)
mkdir -p -- "$2"
work_dir=$(cd -- "$2" && pwd)
workers=${3:-4}
rio_bin=${RIO_BIN:-rio}
partition_bin=${PARTITION_BIN:-partition}
prefix=Bhuvan_CartoDEM_v3r1_TerrainRGB
export GDAL_CACHEMAX=128
export GDAL_NUM_THREADS=1
export TMPDIR="$work_dir/tmp"
mkdir -p "$work_dir/input" "$work_dir/mbtiles" "$work_dir/pmtiles" "$work_dir/logs" "$TMPDIR"

# These are virtual rasters; no extra full-resolution DEM is materialized.
gdalbuildvrt -overwrite -srcnodata -255 -vrtnodata -255 \
    "$work_dir/input/geographic.vrt" "$cog_dir"/*.tif
gdalwarp -overwrite -of VRT -r cubicspline -s_srs EPSG:4326 -t_srs EPSG:3857 \
    -srcnodata -255 -dstnodata 0 \
    "$work_dir/input/geographic.vrt" "$work_dir/input/mercator.vrt"

sources=()
for zoom in {0..12}; do
    round_bits=$((16 - (zoom > 5 ? zoom : 5)))
    output="$work_dir/mbtiles/z${zoom}.mbtiles"
    if [[ ! -f "$output" ]]; then
        partial="$work_dir/mbtiles/z${zoom}.partial.mbtiles"
        if [[ -e "$partial" ]]; then
            echo "Incomplete output exists: $partial; inspect it before restarting." >&2
            exit 1
        fi
        "$rio_bin" rgbify -b -10000 -i 0.1 --encoding mapbox \
            --min-z "$zoom" --max-z "$zoom" --round-digits "$round_bits" \
            -j "$workers" --batch-size 16 --format webp --resampling nearest \
            "$work_dir/input/mercator.vrt" "$partial" \
            >"$work_dir/logs/z${zoom}.log" 2>&1
        python3 "$script_dir/check_mbtiles.py" "$partial"
        mv -- "$partial" "$output"
    else
        python3 "$script_dir/check_mbtiles.py" "$output"
    fi
    sources+=(--from-source "$output")
done

# The existing tool chooses geographic partitions. All thirteen MBTiles sources
# are supplied directly, without exporting/reimporting a tile directory.
"$partition_bin" "${sources[@]}" \
    --to-pmtiles "$work_dir/pmtiles/$prefix.pmtiles" \
    --size-limit 2000000000 --no-cache \
    >"$work_dir/logs/partition.log" 2>&1
