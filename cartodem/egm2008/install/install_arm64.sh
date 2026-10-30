#!/usr/bin/env bash
# Install the pinned processing and partitioning environments on Linux ARM64.
set -euo pipefail

usage() {
    echo 'Usage: bash install/install_arm64.sh TOOLS_DIRECTORY [DEM_WORKSPACE]'
    echo 'Linux ARM64, glibc >= 2.28.'
    echo 'Installs uv, Python and the processing/partitioning dependencies.'
    echo 'Supply DEM_WORKSPACE to install WhiteboxTools for river processing.'
}
if [[ ${1:-} == --help || ${1:-} == -h ]]; then usage; exit 0; fi
if [[ $# -lt 1 || $# -gt 2 ]]; then usage >&2; exit 2; fi
if [[ $(uname -s) != Linux || $(uname -m) != aarch64 ]]; then
    echo 'This installer targets Linux ARM64 (aarch64).' >&2; exit 2
fi
pipeline_code=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
for pipeline_command in curl tar gzip sha256sum; do
    if ! command -v "$pipeline_command" >/dev/null; then
        echo 'Install prerequisites: sudo apt-get update && sudo apt-get install -y ca-certificates curl tar gzip coreutils' >&2
        exit 2
    fi
done
mkdir -p -- "$1"
pipeline_tools=$(cd -- "$1" && pwd)
mkdir -p "$pipeline_tools/bin"
pipeline_scratch=$(mktemp -d "$pipeline_tools/.install.XXXXXXXX")
trap 'rm -rf -- "$pipeline_scratch"' EXIT

# Digest checked against the publisher's release asset metadata, 2026-10-30.
pipeline_uv_version=0.13.6
pipeline_uv_sha=229e194255c7958d86d56e8790d969a6f74736592133f0dfd0076f754c642c1d
pipeline_python_version=3.14.3
curl --fail --silent --show-error --location --retry 3 --connect-timeout 30 --max-time 300 \
    --output "$pipeline_scratch/uv.tar.gz" \
    "https://github.com/astral-sh/uv/releases/download/$pipeline_uv_version/uv-aarch64-unknown-linux-gnu.tar.gz"
printf '%s  %s\n' "$pipeline_uv_sha" "$pipeline_scratch/uv.tar.gz" | sha256sum --check --status
mkdir "$pipeline_scratch/uv"
tar -xzf "$pipeline_scratch/uv.tar.gz" --strip-components=1 -C "$pipeline_scratch/uv"
install -m 755 "$pipeline_scratch/uv/uv" "$pipeline_tools/bin/uv"

# Keep the interpreter, wheel cache and virtual environments on the chosen disk.
export UV_PYTHON_INSTALL_DIR="$pipeline_tools/python"
export UV_CACHE_DIR="$pipeline_tools/cache"
export UV_PYTHON_DOWNLOADS=automatic
"$pipeline_tools/bin/uv" python install "$pipeline_python_version"
for pipeline_env in processing-venv mosaic-venv; do
    if [[ -d "$pipeline_tools/$pipeline_env" ]]; then
        "$pipeline_tools/$pipeline_env/bin/python" -c 'import sys; assert sys.version_info[:3] == (3,14,3), "Use a new installation directory for a different Python version"'
    else
        "$pipeline_tools/bin/uv" venv --python "$pipeline_python_version" "$pipeline_tools/$pipeline_env"
    fi
done
"$pipeline_tools/bin/uv" pip sync --python "$pipeline_tools/processing-venv/bin/python" \
    --only-binary :all: "$pipeline_code/requirements.txt"
(
    cd -- "$pipeline_code"
    "$pipeline_tools/processing-venv/bin/python" -m validation.check_tiles_environment \
        --output "$pipeline_tools/software.json"
)
# pmtiles-mosaic requires a different PMTiles library version.
# Its pinned pmtiles 3.4.1 dependency is distributed as a Python source package.
"$pipeline_tools/bin/uv" pip sync --python "$pipeline_tools/mosaic-venv/bin/python" \
    "$pipeline_code/requirements-mosaic.txt"
"$pipeline_tools/mosaic-venv/bin/partition" --help >/dev/null

if [[ $# == 2 ]]; then
    "$pipeline_tools/processing-venv/bin/python" "$pipeline_code/install/install_whitebox_binary.py" \
        --workdir "$2"
fi
printf '\nProcessing Python: %s\n' "$pipeline_tools/processing-venv/bin/python"
printf 'PMTiles partitioner: %s\n' "$pipeline_tools/mosaic-venv/bin/partition"
echo 'See README.md for downloads, river fitting, COG export and TerrainRGB packaging.'
