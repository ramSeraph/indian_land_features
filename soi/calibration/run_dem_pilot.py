# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "rasterio>=1.4", "scipy>=1.14", "shapely>=2.1", "matplotlib>=3.9"]
# ///
"""Run a limited local SOI DEM pilot and compare every successful sheet to GLO-30."""

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compare_dem import compare
from validate_dem_sheet import file_hash, validate_sheet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sheets", nargs="+", required=True)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/contour_dem/production-pilot")
    )
    parser.add_argument(
        "--comparisons-dir",
        type=Path,
        default=Path("data/comparisons/production-pilot"),
    )
    parser.add_argument(
        "--cache-dir", type=Path, default=Path("data/contours/production-pilot")
    )
    parser.add_argument(
        "--parquet",
        type=Path,
        default=Path("data/contours/soi-raw/SOI_Contours_raw.parquet.meta.json"),
    )
    parser.add_argument(
        "--iomaps", default="../../indianopenmaps/python/.venv/bin/iomaps"
    )
    parser.add_argument("--min-free-gib", type=float, default=3)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    index = root / "index_50k.geojson"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    status_path = args.output_dir / "pilot_status.json"
    status = json.loads(status_path.read_text()) if status_path.exists() else {}
    for sheet in args.sheets:
        if not sheet or Path(sheet).name != sheet or sheet in (".", ".."):
            raise ValueError("Invalid sheet ID")
        if shutil.disk_usage(args.output_dir).free < args.min_free_gib * 1024**3:
            raise RuntimeError(
                "Disk reserve reached; pilot can be resumed after freeing space"
            )
        started = time.monotonic()
        target = args.output_dir / sheet
        comparison = args.comparisons_dir / sheet
        status[sheet] = {"status": "running"}
        status_path.write_text(json.dumps(status, indent=2) + "\n")
        print(f"PILOT {sheet}", flush=True)
        try:
            if not (target / "report.json").exists():
                command = [
                    sys.executable,
                    str(root / "extract_sheet_dem.py"),
                    "--sheet",
                    sheet,
                    "--parquet",
                    str(args.parquet),
                    "--iomaps",
                    args.iomaps,
                    "--output-dir",
                    str(target),
                    "--cache-dir",
                    str(args.cache_dir),
                    "--no-water",
                    "--reuse-extracts",
                ]
                with (args.output_dir / f"{sheet}.log").open("a") as log:
                    subprocess.run(
                        command, check=True, stdout=log, stderr=subprocess.STDOUT
                    )
            checks = validate_sheet(
                target, index, root / "data/contours" / sheet / "contours.geojsonl"
            )
            # The earlier Odisha reference has exactly this native sampling grid.
            # reference_window still verifies its bounds and checksum before reuse.
            if sheet == "73D_15" and not comparison.exists():
                old = root / "data/comparisons/73D_15-glo30-native-grid"
                comparison.mkdir(parents=True)
                for name in ("reference.json", "glo30_native_window.tif"):
                    shutil.copy2(old / name, comparison / name)
            if (comparison / "report.json").exists():
                report = json.loads((comparison / "report.json").read_text())
                assert report["dem_sha256"] == file_hash(target / "dem.tif"), (
                    "Stale GLO-30 comparison"
                )
            else:
                report = compare(
                    target / "dem.tif",
                    comparison,
                    target / "spots.geojsonl",
                    target / "constraint_distance.tif",
                )
            status[sheet] = {
                "status": "compared",
                "sanity_checks": checks,
                "glo30": report["all_common_pixels"],
                "slope_groups": report["slope_groups"],
                "source_points_vs_glo30": report.get("source_points"),
                "glo30_report": str(comparison / "report.json"),
                "elapsed_seconds": time.monotonic() - started,
                "output_bytes": sum(
                    p.stat().st_size for p in target.rglob("*") if p.is_file()
                ),
            }
        except Exception as exc:  # noqa: BLE001 -- record each failed pilot explicitly and continue
            status[sheet] = {
                "status": "failed",
                "error_type": type(exc).__name__,
                "error": str(exc),
                "elapsed_seconds": time.monotonic() - started,
            }
            print(f"FAILED {sheet}: {exc}", flush=True)
        status_path.write_text(json.dumps(status, indent=2, allow_nan=False) + "\n")
        print(f"PILOT RESULT {sheet}: {status[sheet]['status']}", flush=True)


if __name__ == "__main__":
    main()
