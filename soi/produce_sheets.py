# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "rasterio>=1.4", "scipy>=1.14", "shapely>=2.1", "gh_release_tools==0.1.14"]
# ///
"""Resumable SOI sheet production with local overlap correction and verified uploads.

Calibration/holdout and GLO-30 comparisons are separate. The production pass
includes all height points. Water conditioning is deferred to separate research.
No toposheet or reference DEM is used to create elevations.
"""

import argparse
import fcntl
import hashlib
import importlib.metadata
import json
import math
import os
import shutil
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import rasterio
from rasterio.shutil import copy as raster_copy
from release_tools.utils import get_release_map

from download_inputs import checksum, download
from overlap_cache import OverlapProduction
from validate_dem_sheet import validate_sheet

ROOT = Path(__file__).resolve().parent
CODE_FILES = [
    "produce_sheets.py",
    "download_inputs.py",
    "extract_sheet_dem.py",
    "contours_to_dem.py",
    "contour_contacts.py",
    "contour_form_lines.py",
    "contour_height_recovery.py",
    "imagico_dem.py",
    "imagico_interpolate.py",
    "imagico_provenance.json",
    "LICENSE.imagico.txt",
    "validate_dem_sheet.py",
    "sheet_overlaps.py",
    "seam_interpolate.py",
    "overlap_cache.py",
]


def atomic_json(path, value):
    staged = path.with_name(path.name + ".tmp")
    staged.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    staged.replace(path)


def digest_object(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def code_hashes():
    return {name: checksum(ROOT / name) for name in CODE_FILES}


def iomaps_identity(executable):
    path = Path(shutil.which(executable) or executable).resolve()
    python = path.parent / "python"
    query = (
        "import hashlib,importlib.metadata,json,pathlib,iomaps; "
        "root=pathlib.Path(iomaps.__file__).parent; "
        "names=['commands/extract.py','commands/decorators.py','core/local_parquet.py','core/spatial_filter.py','core/helpers.py']; "
        "print(json.dumps({'version':importlib.metadata.version('iomaps'),'files':{n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in names}},sort_keys=True))"
    )
    return json.loads(subprocess.check_output([str(python), "-c", query], text=True))


def package_sheet(sheet, job, destination, outcome):
    """Export only the elevation COG; retain requested height points separately."""
    destination.mkdir(parents=True, exist_ok=True)
    dem_dir = job / "dem"
    basename = f"SOI_DEM_30m_{sheet}"
    if outcome["status"] == "built":
        cog = destination / f"{basename}.tif"
        raster_copy(
            dem_dir / "dem.tif",
            cog,
            driver="COG",
            compress="DEFLATE",
            predictor="FLOATING_POINT",
            blocksize=512,
            overview_resampling="AVERAGE",
            num_threads="1",
        )
        with rasterio.open(dem_dir / "dem.tif") as source, rasterio.open(cog) as target:
            if (source.crs, source.transform, source.shape, source.nodata) != (
                target.crs,
                target.transform,
                target.shape,
                target.nodata,
            ):
                raise ValueError("COG export changed the grid")
            if not np.array_equal(source.read(1), target.read(1)):
                raise ValueError("COG export changed elevations")
            if target.tags(ns="IMAGE_STRUCTURE").get("LAYOUT") != "COG":
                raise ValueError("Expected cloud-optimized GeoTIFF")
    return [
        {"name": p.name, "bytes": p.stat().st_size, "sha256": checksum(p)}
        for p in sorted(destination.iterdir())
        if p.is_file()
    ]


def retain_height_points(sheet, dem_dir, destination):
    source = dem_dir / "spots.geojsonl"
    if source.exists():
        destination.mkdir(parents=True, exist_ok=True)
        staged = destination / f"{sheet}.geojsonl.tmp"
        shutil.copy2(source, staged)
        staged.replace(destination / f"{sheet}.geojsonl")


def upload_and_verify(folder, assets, repo, release, receipt):
    # gh-release-tools owns release overflow and upload batching. It skips by
    # name, so independently verify GitHub's digests before releasing local space.
    executable = shutil.which("upload-to-release")
    if executable is None:
        raise ValueError(
            "Run with the locked uv environment containing gh-release-tools"
        )
    subprocess.run(
        [
            executable,
            "--repo",
            repo,
            "--release",
            release,
            "--folder",
            str(folder),
            "--extension",
            ".tif",
            "--create-extra-releases",
            "--batch-size",
            "1",
        ],
        check=True,
    )
    wanted = {a["name"]: a for a in assets}
    verified = []
    for tag in get_release_map(release, repo=repo).values():
        remote = json.loads(
            subprocess.check_output(
                ["gh", "api", f"repos/{repo}/releases/tags/{tag}"], text=True
            )
        )
        pages = json.loads(
            subprocess.check_output(
                [
                    "gh",
                    "api",
                    "--paginate",
                    "--slurp",
                    f"repos/{repo}/releases/{remote['id']}/assets?per_page=100",
                ],
                text=True,
            )
        )
        for asset in (a for page in pages for a in page):
            name = asset["name"]
            if name not in wanted:
                continue
            expected = wanted[name]
            if (
                asset["size"] != expected["bytes"]
                or asset.get("digest") != "sha256:" + expected["sha256"]
            ):
                raise ValueError(
                    f"Remote asset differs from local output: {tag}/{name}"
                )
            verified.append(
                {
                    **expected,
                    "asset_id": asset["id"],
                    "release": tag,
                    "url": asset["browser_download_url"],
                }
            )
            del wanted[name]
    if wanted:
        raise ValueError(f"Upload missing assets: {sorted(wanted)}")
    result = {
        "status": "verified",
        "verified_at": datetime.now(UTC).isoformat(),
        "repo": repo,
        "assets": verified,
    }
    atomic_json(receipt, result)
    return result


def build_independent(sheet, job, args, state, state_path):
    started = time.monotonic()
    command = [
        sys.executable,
        str(ROOT / "extract_sheet_dem.py"),
        "--sheet",
        sheet,
        "--parquet",
        str(args.parquet.resolve()),
        "--iomaps",
        args.iomaps,
        "--output-dir",
        str(job / "dem"),
        "--cache-dir",
        str(job / "source"),
        "--no-water",
        "--skip-holdout",
    ]
    try:
        with (
            (job / "build.log").open("w") as log,
            subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT) as child,
        ):
            state["child_pid"] = child.pid
            atomic_json(state_path, state)
            try:
                code = child.wait(timeout=args.sheet_timeout)
            except (subprocess.TimeoutExpired, KeyboardInterrupt):
                child.kill()
                child.wait()
                raise
            finally:
                state.pop("child_pid", None)
                atomic_json(state_path, state)
            if code:
                raise subprocess.CalledProcessError(code, command)
        checks = validate_sheet(job / "dem", ROOT / "index_50k.geojson")
        result = {"status": "built", "sheet": sheet, "checks": checks}
    except (
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
        ValueError,
        AssertionError,
        OSError,
    ) as exc:
        result = {
            "status": "failed",
            "sheet": sheet,
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
    return {**result, "initial_elapsed_seconds": time.monotonic() - started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--sheets-file",
        type=Path,
        default=ROOT / "data/contours/soi-raw/sheet-inventory/sheets_interior.txt",
    )
    parser.add_argument(
        "--parquet",
        type=Path,
        default=ROOT / "data/contours/soi-raw/SOI_Contours_raw.parquet.meta.json",
    )
    parser.add_argument(
        "--inputs-manifest", type=Path, default=ROOT / "manifests/soi-contours.json"
    )
    parser.add_argument(
        "--iomaps", default=str(ROOT / "../../indianopenmaps/python/.venv/bin/iomaps")
    )
    parser.add_argument("--priority", nargs="*", default=[])
    parser.add_argument(
        "--limit",
        type=int,
        help="Process this many additional sheets, then stop cleanly",
    )
    parser.add_argument("--min-free-gib", type=float, default=3)
    parser.add_argument("--sheet-timeout", type=float, default=1800)
    parser.add_argument(
        "--no-overlap",
        action="store_true",
        help="Produce independent first-pass sheets only",
    )
    parser.add_argument("--overlap-width-m", type=float, default=2000)
    parser.add_argument("--overlap-context-m", type=float, default=3600)
    parser.add_argument("--repo", default="ramSeraph/indian_land_features")
    parser.add_argument(
        "--release",
        help="Upload with gh-release-tools to an existing validated release",
    )
    parser.add_argument(
        "--prune-uploaded",
        action="store_true",
        help="Remove this run's local packages only after remote size and SHA256 verification",
    )
    args = parser.parse_args()
    if args.prune_uploaded and not args.release:
        parser.error("--prune-uploaded requires --release")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if args.overlap_width_m <= 0 or args.overlap_context_m < 3600:
        parser.error(
            "Overlap requires a positive width and at least 3,600 m of context"
        )
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    lock = (out / ".lock").open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        parser.error("Another producer is using this output directory")
    sheets = args.sheets_file.read_text().splitlines()
    if len(set(sheets)) != len(sheets) or any(
        not s or s in (".", "..") or Path(s).name != s for s in sheets
    ):
        parser.error("Sheet list must contain distinct safe IDs")
    if any(s not in sheets for s in args.priority):
        parser.error("Priority sheets must belong to the selected sheet list")
    # Spatial order keeps adjacent first passes in the bounded cache longer.
    indexed = {
        f["properties"]["id"]: f["geometry"]["coordinates"][0]
        for f in json.loads((ROOT / "index_50k.geojson").read_text())["features"]
    }
    spatial = sorted(
        sheets,
        key=lambda s: (
            -math.floor(max(p[1] for p in indexed[s]) * 2 - 1e-7),
            math.floor(min(p[0] for p in indexed[s]) * 2 + 1e-7),
            -max(p[1] for p in indexed[s]),
            min(p[0] for p in indexed[s]),
            s,
        ),
    )
    ordered = list(dict.fromkeys(args.priority + spatial))
    # The run manifest retains the verified input identities.
    receipt = download(
        args.inputs_manifest, args.parquet.resolve().parent, args.min_free_gib
    )
    config = {
        "schema_version": 3,
        "output_products": ["elevation_cog"],
        "resolution_m": 30,
        "interpolation": "imagico",
        "sheet_grid_crs": "EPSG:4326",
        "validate_points": False,
        "water_conditioning": "deferred_to_separate_research",
        "overlap": {
            "enabled": not args.no_overlap,
            "width_per_side_m": args.overlap_width_m,
            "context_m": args.overlap_context_m,
            "corner_shape": "intersection_of_incident_edge_strips",
            "source_cache_sheets": 16,
            "patch_cache_entries": 64,
        },
        "index_sha256": checksum(ROOT / "index_50k.geojson"),
        "input_manifest_sha256": checksum(args.inputs_manifest),
        "sheets": sheets,
        "code": code_hashes(),
        "packages": {
            p: importlib.metadata.version(p)
            for p in (
                "numpy",
                "scipy",
                "rasterio",
                "pyproj",
                "shapely",
                "gh_release_tools",
            )
        },
        "python": sys.version,
        "gdal": rasterio.__gdal_version__,
        "iomaps": iomaps_identity(args.iomaps),
    }
    fingerprint = digest_object(config)
    plan_path, state_path = out / "run.json", out / "state.json"
    if plan_path.exists():
        previous = json.loads(plan_path.read_text())
        if previous["fingerprint"] != fingerprint:
            raise ValueError(
                "Code, input, parameters or environment changed; choose a new output directory"
            )
    else:
        atomic_json(
            plan_path,
            {
                "fingerprint": fingerprint,
                "configuration": config,
                "input_verification": receipt,
            },
        )
    state = (
        json.loads(state_path.read_text())
        if state_path.exists()
        else {"fingerprint": fingerprint, "sheets": {}}
    )
    if state["fingerprint"] != fingerprint:
        raise ValueError("Checkpoint belongs to a different run")
    if state.get("child_pid"):
        try:
            os.kill(state["child_pid"], 0)
        except ProcessLookupError:
            state.pop("child_pid", None)
        else:
            raise ValueError("A child from the previous producer may still be running")
    for directory in ("jobs", "packages", "receipts", "reports", "logs"):
        (out / directory).mkdir(exist_ok=True)
    processed = 0
    overlap = (
        None
        if args.no_overlap
        else OverlapProduction(
            out,
            ROOT / "index_50k.geojson",
            sheets,
            fingerprint,
            lambda sheet, job: build_independent(sheet, job, args, state, state_path),
            width_m=args.overlap_width_m,
            context_m=args.overlap_context_m,
            reserve_gib=args.min_free_gib,
        )
    )
    try:
        for sheet in ordered:
            old = state["sheets"].get(sheet, {})
            if old.get("status") in {"uploaded", "failed"}:
                continue
            if old.get("status") == "ready" and not args.release:
                continue
            if (out / "STOP").exists():
                state["status"] = "stopped_by_request"
                break
            if (
                overlap is not None
                and shutil.disk_usage(out).free < args.min_free_gib * 1024**3
            ):
                try:
                    overlap.space_check()
                except OSError:
                    state["status"] = "waiting_for_disk_space"
                    break
            if shutil.disk_usage(out).free < args.min_free_gib * 1024**3:
                state["status"] = "waiting_for_disk_space"
                break
            if code_hashes() != config["code"]:
                raise ValueError(
                    "Code changed during processing; stopping before the next sheet"
                )
            for asset in receipt["assets"]:
                stat = Path(asset["path"]).stat()
                if (stat.st_size, stat.st_mtime_ns) != (
                    asset["bytes"],
                    asset["mtime_ns"],
                ):
                    raise ValueError("Source input changed during processing")
            print(f"SHEET {sheet}", flush=True)
            state.update(status="running", current_sheet=sheet, pid=os.getpid())
            atomic_json(state_path, state)
            destination = out / "packages" / sheet
            if old.get("status") == "ready":
                assets = old["assets"]
                for asset in assets:
                    p = destination / asset["name"]
                    if (
                        not p.is_file()
                        or p.stat().st_size != asset["bytes"]
                        or checksum(p) != asset["sha256"]
                    ):
                        raise ValueError(f"Saved package changed: {p}")
                outcome = old["outcome"]
            else:
                job = out / "jobs" / sheet
                if job.exists() or destination.exists():
                    child_pid = state.get("child_pid")
                    if child_pid:
                        try:
                            os.kill(child_pid, 0)
                        except ProcessLookupError:
                            pass
                        else:
                            raise ValueError(
                                f"Previous sheet process {child_pid} may still be running; waiting for it before resuming"
                            )
                    interrupted = out / "interrupted" / f"{sheet}-{time.time_ns()}"
                    interrupted.mkdir(parents=True)
                    if job.exists():
                        job.replace(interrupted / "job")
                    if destination.exists():
                        destination.replace(interrupted / "packages")
                job.mkdir()
                started = time.monotonic()
                if overlap is None:
                    outcome = build_independent(sheet, job, args, state, state_path)
                else:
                    outcome = overlap.prepare(sheet, job)
                    if outcome["status"] == "built":
                        outcome = {
                            **outcome,
                            "checks": validate_sheet(
                                job / "dem", ROOT / "index_50k.geojson"
                            ),
                        }
                outcome.update(
                    elapsed_seconds=time.monotonic() - started,
                    run_fingerprint=fingerprint,
                )
                if code_hashes() != config["code"]:
                    raise ValueError(
                        "Code changed while building this sheet; refusing to package it"
                    )
                for asset in receipt["assets"]:
                    stat = Path(asset["path"]).stat()
                    if (stat.st_size, stat.st_mtime_ns) != (
                        asset["bytes"],
                        asset["mtime_ns"],
                    ):
                        raise ValueError(
                            "Source input changed while building this sheet"
                        )
                assets = package_sheet(sheet, job, destination, outcome)
                retain_height_points(sheet, job / "dem", ROOT / "data/height_points")
                atomic_json(out / "reports" / f"{sheet}.json", outcome)
                shutil.copy2(job / "build.log", out / "logs" / f"{sheet}.log")
                # Inputs remain reproducible from the pinned source data. Keep
                # the requested point dataset and compact run reports, then
                # release this job's temporary extraction/interpolation files.
                shutil.rmtree(job)
                state["sheets"][sheet] = {
                    "status": "ready" if assets else "failed",
                    "outcome": outcome,
                    "assets": assets,
                }
                atomic_json(state_path, state)
            if args.release and assets:
                published = upload_and_verify(
                    destination,
                    assets,
                    args.repo,
                    args.release,
                    out / "receipts" / f"{sheet}.json",
                )
                state["sheets"][sheet].update(status="uploaded", publication=published)
                atomic_json(state_path, state)
                if args.prune_uploaded:
                    shutil.rmtree(destination)
            if overlap is not None:
                overlap.trim()
            processed += 1
            print(
                f"RESULT {sheet}: {outcome['status']}; {state['sheets'][sheet]['status']}",
                flush=True,
            )
            if args.limit is not None and processed >= args.limit:
                state["status"] = "limit_reached"
                break
        else:
            state["status"] = "finished"
    except BaseException as exc:
        state.update(
            status="interrupted" if isinstance(exc, KeyboardInterrupt) else "error",
            error_type=type(exc).__name__,
            error=str(exc),
        )
        raise
    finally:
        state["updated_at"] = datetime.now(UTC).isoformat()
        state["counts"] = dict(
            __import__("collections").Counter(
                r["outcome"]["status"] for r in state["sheets"].values()
            )
        )
        atomic_json(state_path, state)
    print(
        json.dumps({"status": state["status"], "counts": state["counts"]}), flush=True
    )


if __name__ == "__main__":
    main()
