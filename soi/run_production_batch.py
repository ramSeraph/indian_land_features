"""Run the resumable unconditioned DEM queue and refresh its release listing."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--sheets-file", type=Path, default=root / "manifests/sheets_interior.txt"
    )
    parser.add_argument("--repo", default="ramSeraph/indian_land_features")
    parser.add_argument("--release", required=True)
    parser.add_argument("--priority", nargs="*", default=[])
    args = parser.parse_args()
    command = [
        "uv",
        "run",
        "--locked",
        str(root / "produce_sheets.py"),
        "--output-dir",
        str(args.output_dir.resolve()),
        "--sheets-file",
        str(args.sheets_file.resolve()),
        "--repo",
        args.repo,
        "--release",
        args.release,
        "--prune-uploaded",
    ]
    if args.priority:
        command.extend(["--priority", *args.priority])
    subprocess.run(command, cwd=root, check=True)
    state = json.loads((args.output_dir / "state.json").read_text())
    if state["status"] != "finished":
        print("Queue stopped without completing: " + state["status"], flush=True)
        return
    subprocess.run(
        [
            "uv",
            "run",
            "--with",
            "gh_release_tools==0.1.14",
            "generate-lists",
            "--repo",
            args.repo,
            "--release",
            args.release,
            "--extension",
            ".tif",
        ],
        cwd=root,
        check=True,
    )
    print(
        "Queue complete; release listing refreshed. " + json.dumps(state["counts"]),
        flush=True,
    )


if __name__ == "__main__":
    main()
