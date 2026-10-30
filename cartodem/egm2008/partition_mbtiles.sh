#!/usr/bin/env bash
# Usage: partition_mbtiles.sh MBTILES_DIRECTORY PMTILES_DIRECTORY [PARTITION_COMMAND]
# The partition command is supplied by pmtiles-mosaic==0.2.0.
set -euo pipefail
if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 MBTILES_DIRECTORY PMTILES_DIRECTORY [PARTITION_COMMAND]" >&2
    exit 2
fi
source_dir=$(cd -- "$1" && pwd)
mkdir -p -- "$2"
output_dir=$(cd -- "$2" && pwd)
partition_command=${3:-partition}
mkdir -p "$output_dir/tmp"
export TMPDIR="$output_dir/tmp"
sources=()
for zoom in {0..12}; do
    source_file="$source_dir/z${zoom}.mbtiles"
    [[ -f "$source_file" ]] || { echo "Missing source: $source_file" >&2; exit 1; }
    sources+=(--from-source "$source_file")
done
exec "$partition_command" "${sources[@]}" \
    --to-pmtiles "$output_dir/Bhuvan_CartoDEM_v3r1_TerrainRGB.pmtiles" \
    --size-limit 2000000000 --no-cache
