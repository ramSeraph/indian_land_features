#!/usr/bin/env bash
# Install the complete pinned DEM, tiling and partitioning environments.
set -euo pipefail
if [[ $# != 2 ]]; then
    echo "Usage: bash $0 TOOLS_DIRECTORY DEM_WORKSPACE" >&2
    exit 2
fi
pipeline_code=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
bash "$pipeline_code/install_tiles_arm64.sh" "$1"
pipeline_tools=$(cd -- "$1" && pwd)
mkdir -p -- "$2"
pipeline_work=$(cd -- "$2" && pwd)
export UV_PYTHON_INSTALL_DIR="$pipeline_tools/python"
export UV_CACHE_DIR="$pipeline_tools/cache"
if [[ ! -d "$pipeline_tools/processing-venv" ]]; then
    "$pipeline_tools/bin/uv" venv --python "$pipeline_tools/venv/bin/python" \
        "$pipeline_tools/processing-venv"
fi
"$pipeline_tools/bin/uv" pip sync --python "$pipeline_tools/processing-venv/bin/python" \
    --only-binary :all: "$pipeline_code/requirements.txt"
"$pipeline_tools/processing-venv/bin/python" "$pipeline_code/install_whitebox_binary.py" \
    --workdir "$pipeline_work"
printf '\nProcessing Python: %s\n' "$pipeline_tools/processing-venv/bin/python"
echo 'See README.md for downloads, river fitting, COG export and TerrainRGB packaging.'
