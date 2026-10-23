"""Run the resumable unconditioned DEM queue and refresh its release listing."""

import argparse
import json
import subprocess
import time
from pathlib import Path


def run_producer(command, root, output_dir):
    """Resume a checked, ready-to-upload sheet after an upload command fails."""
    state_path = output_dir / "state.json"
    retry_sheet, retries = None, 0
    delays = (15, 30, 60, 120)
    while True:
        before = state_path.stat().st_mtime_ns if state_path.exists() else None
        try:
            subprocess.run(command, cwd=root, check=True)
            return
        except subprocess.CalledProcessError:
            # Startup/configuration failures leave the old checkpoint untouched.
            # Do not misread an earlier upload error as the current failure.
            if not state_path.exists() or state_path.stat().st_mtime_ns == before:
                raise
            state = json.loads(state_path.read_text())
            sheet = state.get("current_sheet")
            record = state.get("sheets", {}).get(sheet, {})
            if (
                state.get("status") != "error"
                or state.get("error_type") != "CalledProcessError"
                or state.get("child_pid")
                or record.get("status") != "ready"
                or record.get("outcome", {}).get("status") != "built"
            ):
                raise
            if sheet != retry_sheet:
                retry_sheet, retries = sheet, 0
            if retries == len(delays):
                raise
            delay = delays[retries]
            retries += 1
            print(
                f"Upload command failed for {sheet}; resuming its saved COG "
                f"in {delay}s (retry {retries}/{len(delays)}).",
                flush=True,
            )
            time.sleep(delay)


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
    run_producer(command, root, args.output_dir.resolve())
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
