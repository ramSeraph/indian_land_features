#!/usr/bin/env bash
# Install the COG-to-TerrainRGB toolchain on Linux ARM64, without rebuilding DEMs.
set -euo pipefail

usage() {
    echo 'Usage: bash install_tiles_arm64.sh /path/to/tiling-tools'
    echo 'Linux ARM64, glibc >= 2.28 (Ubuntu 22.04/24.04 or Debian 12 recommended).'
    echo 'Installs pinned uv, Python, tiling wheels and GitHub CLI into this directory.'
}
if [[ ${1:-} == --help || ${1:-} == -h ]]; then usage; exit 0; fi
if [[ $# != 1 ]]; then usage >&2; exit 2; fi
if [[ $(uname -s) != Linux || $(uname -m) != aarch64 ]]; then
    echo 'This installer targets Linux ARM64 (aarch64).' >&2; exit 2
fi
tile_code_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
for tile_command in curl tar gzip sha256sum; do
    if ! command -v "$tile_command" >/dev/null; then
        echo 'Install prerequisites: sudo apt-get update && sudo apt-get install -y ca-certificates curl tar gzip coreutils' >&2
        exit 2
    fi
done
mkdir -p -- "$1"
tile_prefix=$(cd -- "$1" && pwd)
mkdir -p "$tile_prefix/bin"
tile_scratch=$(mktemp -d "$tile_prefix/.install.XXXXXXXX")
trap 'rm -rf -- "$tile_scratch"' EXIT

# Digests are from the publishers release asset metadata, checked 2026-10-30.
tile_uv_version=0.13.6
tile_uv_sha=229e194255c7958d86d56e8790d969a6f74736592133f0dfd0076f754c642c1d
tile_gh_version=2.103.2
tile_gh_sha=4c3e99a50bd8b181a3fe9d69b25dd978bac768148d88f01507c6b4d1ef61a443
tile_python_version=3.14.3

download_verified() {
    curl --fail --location --retry 3 --output "$tile_scratch/$2" "$1"
    printf '%s  %s\n' "$3" "$tile_scratch/$2" | sha256sum --check --status
}
download_verified "https://github.com/astral-sh/uv/releases/download/$tile_uv_version/uv-aarch64-unknown-linux-gnu.tar.gz" uv.tar.gz "$tile_uv_sha"
mkdir "$tile_scratch/uv"
tar -xzf "$tile_scratch/uv.tar.gz" --strip-components=1 -C "$tile_scratch/uv"
install -m 755 "$tile_scratch/uv/uv" "$tile_prefix/bin/uv"

download_verified "https://github.com/cli/cli/releases/download/v$tile_gh_version/gh_${tile_gh_version}_linux_arm64.tar.gz" gh.tar.gz "$tile_gh_sha"
mkdir "$tile_scratch/gh"
tar -xzf "$tile_scratch/gh.tar.gz" --strip-components=1 -C "$tile_scratch/gh"
install -m 755 "$tile_scratch/gh/bin/gh" "$tile_prefix/bin/gh"

# Keep the interpreter, wheel cache and virtual environment on the chosen disk.
export UV_PYTHON_INSTALL_DIR="$tile_prefix/python"
export UV_CACHE_DIR="$tile_prefix/cache"
export UV_PYTHON_DOWNLOADS=automatic
"$tile_prefix/bin/uv" python install "$tile_python_version"
if [[ -d "$tile_prefix/venv" ]]; then
    "$tile_prefix/venv/bin/python" -c 'import sys; assert sys.version_info[:3] == (3,14,3), "Use a new installation directory for a different Python version"'
else
    "$tile_prefix/bin/uv" venv --python "$tile_python_version" "$tile_prefix/venv"
fi
"$tile_prefix/bin/uv" pip sync --python "$tile_prefix/venv/bin/python" \
    --only-binary :all: "$tile_code_dir/requirements-tiles.txt"
(
    cd -- "$tile_code_dir"
    "$tile_prefix/venv/bin/python" -m validation.check_tiles_environment \
        --output "$tile_prefix/software.json"
)
# pmtiles-mosaic pins a different PMTiles library version, so keep its environment separate.
if [[ ! -d "$tile_prefix/mosaic-venv" ]]; then
    "$tile_prefix/bin/uv" venv --python "$tile_python_version" "$tile_prefix/mosaic-venv"
fi
"$tile_prefix/bin/uv" pip sync --python "$tile_prefix/mosaic-venv/bin/python" \
    --only-binary :all: "$tile_code_dir/requirements-mosaic.txt"
"$tile_prefix/mosaic-venv/bin/partition" --help >/dev/null
"$tile_prefix/bin/uv" --version
"$tile_prefix/bin/gh" --version
printf '\nReady. Python: %s\nGitHub CLI: %s\n' "$tile_prefix/venv/bin/python" "$tile_prefix/bin/gh"
printf 'PMTiles partitioner: %s\n' "$tile_prefix/mosaic-venv/bin/partition"
echo 'Authenticate GitHub CLI separately with: gh auth login --hostname github.com --git-protocol https --web'
echo 'See TILING_SETUP.md for rebuilding from existing COGs and uploading release assets.'
