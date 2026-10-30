#!/usr/bin/env bash
# Encode one final-resolution numeric elevation tile using GDAL command-line tools.
# Usage: bash gdal_encode_tile.sh heights.tif output.webp ROUND_BITS
# Input: one elevation band in metres, NaN nodata, already resampled for this zoom.
set -euo pipefail

if [[ $# != 3 ]]; then
    echo 'Usage: bash gdal_encode_tile.sh heights.tif output.webp ROUND_BITS' >&2
    exit 2
fi
source_path=$1
output_path=$2
round_bits=$3
if [[ ! $round_bits =~ ^[0-9]+$ ]] || (( 10#$round_bits > 23 )); then
    echo 'ROUND_BITS must be an integer from 0 to 23' >&2
    exit 2
fi
if [[ -e $output_path ]]; then
    echo "Refusing to overwrite $output_path" >&2
    exit 2
fi
command -v gdal_calc.py >/dev/null
command -v gdal_translate >/dev/null
output_dir=$(dirname "$output_path")
mkdir -p "$output_dir"
scratch=$(mktemp -d "$output_dir/.gdal-encode.XXXXXX")
trap 'rm -rf -- "$scratch"' EXIT

step=$((1 << 10#$round_bits))
# Replace invalid values only while calculating their unused RGB code. Each
# output band below explicitly restores the exact nodata RGBA placeholder.
heights='nan_to_num(A.astype(float64),nan=0,posinf=0,neginf=0)'
code="(rint((($heights+10000)/0.1)/$step)*$step).astype(uint32)"
gdal_calc.py -A "$source_path" --hideNoData --type=Byte --NoDataValue=none \
    --calc="where(isfinite(A),($code>>16)&255,1)" \
    --calc="where(isfinite(A),($code>>8)&255,134)" \
    --calc="where(isfinite(A),$code&255,160)" \
    --calc='full_like(A,255)' \
    --co TILED=YES --co COMPRESS=DEFLATE \
    --outfile="$scratch/rgba.tif" --quiet
gdal_translate -q -of WEBP -co LOSSLESS=YES \
    -colorinterp red,green,blue,alpha "$scratch/rgba.tif" "$scratch/tile.webp"
# Copy with an exclusive destination: do not overwrite a file that appeared
# during encoding. The temporary directory is removed by the EXIT trap.
if ! (set -o noclobber; cat "$scratch/tile.webp" > "$output_path"); then
    echo "Could not create $output_path without overwriting it" >&2
    exit 1
fi
