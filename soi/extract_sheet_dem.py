#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "rasterio>=1.4", "scipy>=1.14", "shapely>=2.1"]
# ///
"""Extract an SOI sheet with iomaps and build its DEM without a scanned toposheet."""

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from urllib.request import urlopen

from contour_form_lines import parse_form_line_interval
from contours_to_dem import build_dem, geometry_hash, get_sheet, read_features

CONTOURS = "/dem/contours/raw/soi/"
WATER = "/dem/water-bodies/cartodem-v3r1/bhuvan/"
INDEX_URL = (
    "https://github.com/ramSeraph/india_topo_maps/releases/download/"
    "soi-ancillary/index_50k.geojson"
)


def ensure_index(index):
    """Fetch the default sheet index if absent; use a supplied index as-is."""
    if index is not None:
        if not index.is_file():
            raise ValueError(f"Sheet index does not exist: {index}")
        return index
    index = Path(__file__).with_name("index_50k.geojson")
    if not index.exists():
        with tempfile.TemporaryDirectory(dir=index.parent) as temporary:
            download = Path(temporary) / index.name
            with urlopen(INDEX_URL, timeout=60) as response:
                download.write_bytes(response.read())
            data = json.loads(download.read_text())
            if data.get("type") != "FeatureCollection" or not data.get("features"):
                raise ValueError("Downloaded sheet index is not a FeatureCollection")
            download.replace(index)
    return index


def validate_extract(path, source):
    expected = (
        {"LineString", "MultiLineString", "Point"}
        if source == CONTOURS
        else {"Polygon", "MultiPolygon"}
    )
    count = 0
    for feature in read_features(path):
        geometry = feature.get("geometry") or {}
        if geometry.get("type") not in expected:
            raise ValueError(f"Unexpected geometry in {path}: {geometry.get('type')}")
        count += 1
    if source == CONTOURS and count == 0:
        raise ValueError("No SOI contours found for this sheet; cannot build a DEM")
    return count


def extract_layer(
    iomaps,
    source,
    sheet_id,
    sheet_hash,
    index,
    target,
    reuse=False,
    http_client="default",
    local_parquets=None,
):
    metadata_path = target.with_suffix(".metadata.json")
    identity = {
        "sheet": sheet_id,
        "source": source,
        "sheet_geometry_sha256": sheet_hash,
    }
    if local_parquets:
        local_files = set()
        for value in local_parquets:
            path = Path(value).resolve()
            local_files.add(path)
            if path.name.endswith(".parquet.meta.json"):
                for name in json.loads(path.read_text())["extents"]:
                    if Path(name).name != name:
                        raise ValueError("Invalid local partition name")
                    local_files.add(path.parent / name)
        identity["local_inputs"] = [
            {
                "path": str(p),
                "bytes": p.stat().st_size,
                "mtime_ns": p.stat().st_mtime_ns,
            }
            for p in sorted(local_files)
        ]
    if target.exists():
        if not reuse:
            raise ValueError(f"Extract already exists: {target}; use --reuse-extracts")
        metadata = (
            json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
        )
        checksum = hashlib.sha256(target.read_bytes()).hexdigest()
        if metadata and (
            any(metadata.get(key) != value for key, value in identity.items())
            or metadata.get("sha256") != checksum
        ):
            raise ValueError(f"Cached extract does not match its metadata: {target}")
        count = validate_extract(target, source)
        print(f"Using {target} ({count} features)", flush=True)
        return {
            **identity,
            "path": str(target),
            "sha256": checksum,
            "features": count,
            "reused": True,
            "origin": metadata.get(
                "origin", "existing_extract_without_provenance_manifest"
            ),
        }
    executable = shutil.which(iomaps)
    if executable is None:
        raise ValueError(f"Cannot find iomaps executable: {iomaps}; use --iomaps PATH")
    target.parent.mkdir(parents=True, exist_ok=True)
    # iomaps/GDAL may leave an incomplete output after an interrupted transfer.
    # Only a successful, validated result becomes a reusable extract.
    with tempfile.TemporaryDirectory(prefix="extract-", dir=target.parent) as temporary:
        staged = Path(temporary) / target.name
        command = [executable, "cli", "extract"]
        if local_parquets:
            for path in local_parquets:
                command.extend(["--parquet", str(Path(path).resolve())])
        else:
            command.extend(["-s", source])
        command.extend(
            [
                "-f",
                str(index),
                "--pick-filter-feature-kv",
                f"id={sheet_id}",
                "-o",
                str(staged),
            ]
        )
        if http_client != "default" and not local_parquets:
            command.extend(["--http-client", http_client])
        subprocess.run(command, check=True)
        count = validate_extract(staged, source)
        checksum = hashlib.sha256(staged.read_bytes()).hexdigest()
        staged.replace(target)
    metadata = {
        **identity,
        "path": str(target),
        "sha256": checksum,
        "features": count,
        "reused": False,
        "origin": "iomaps_local_parquet" if local_parquets else "iomaps",
        "http_client": http_client,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def run_sheet(
    sheet_id,
    index,
    output_dir,
    cache_dir,
    iomaps="iomaps",
    reuse=False,
    water=True,
    resolution=30.0,
    spot_radius=13.0,
    triangle_side=70.0,
    http_client="default",
    edge_sliver_m=2.0,
    form_line_interval="auto",
    contour_contact_radius_m=None,
    height_recovery="auto",
    height_repair_max_length_m=10.0,
    height_repair_max_fraction=0.10,
    interpolation="imagico",
    local_parquets=None,
    validate_points=True,
):
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("Output directory is not empty; choose a new run directory")
    sheet = get_sheet(index, sheet_id)
    sheet_hash = geometry_hash(sheet)
    # A sheet ID is also used as a directory name; do not allow path traversal.
    if not sheet_id or sheet_id in {".", ".."} or Path(sheet_id).name != sheet_id:
        raise ValueError("Sheet ID must be a single directory name")
    cache = cache_dir / sheet_id
    sources = [(CONTOURS, cache / "contours.geojsonl")]
    if water:
        sources.append((WATER, cache / "water.geojsonl"))
    extraction = [
        extract_layer(
            iomaps,
            source,
            sheet_id,
            sheet_hash,
            index,
            path,
            reuse,
            http_client,
            local_parquets=local_parquets if source == CONTOURS else None,
        )
        for source, path in sources
    ]
    report = build_dem(
        [sources[0][1]],
        index,
        sheet_id,
        output_dir,
        resolution=resolution,
        water=sources[1][1] if water else None,
        spot_radius=spot_radius,
        triangle_side=triangle_side,
        edge_sliver_m=edge_sliver_m,
        form_line_interval=form_line_interval,
        contour_contact_radius_m=contour_contact_radius_m,
        height_recovery=height_recovery,
        height_repair_max_length_m=height_repair_max_length_m,
        height_repair_max_fraction=height_repair_max_fraction,
        interpolation=interpolation,
        validate_points=validate_points,
    )
    (output_dir / "extraction.json").write_text(
        json.dumps({"sheet": sheet_id, "layers": extraction}, indent=2) + "\n"
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sheet", required=True)
    parser.add_argument(
        "--index", type=Path, help="Default: local index, downloaded if absent"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, default=Path("data/contours"))
    parser.add_argument(
        "--iomaps", default="iomaps", help="Path to the patched iomaps executable"
    )
    parser.add_argument(
        "--reuse-extracts",
        action="store_true",
        help="Reuse existing local extracts after validation",
    )
    parser.add_argument(
        "--no-water",
        action="store_true",
        help="Omit water-mask extraction and smoothing",
    )
    parser.add_argument("--resolution", type=float, default=30.0)
    parser.add_argument(
        "--skip-holdout",
        action="store_true",
        help="Omit calibration-only point holdout; final interpolation still includes all points",
    )
    parser.add_argument(
        "--parquet",
        action="append",
        type=Path,
        help="Local contour Parquet or partition manifest passed to iomaps; repeat for multiple files",
    )
    parser.add_argument(
        "--interpolation", choices=["imagico", "linear"], default="imagico"
    )
    parser.add_argument("--spot-radius", type=float, default=13.0)
    parser.add_argument("--triangle-side", type=float, default=70.0)
    parser.add_argument("--edge-sliver-m", type=float, default=2.0)
    parser.add_argument(
        "--contour-contact-radius-m",
        type=float,
        help="Relax local continuing-contour contacts; default: half grid spacing; 0 disables",
    )
    parser.add_argument(
        "--form-line-interval",
        type=parse_form_line_interval,
        default="auto",
        help="auto (default): infer one interval per sheet; none: disable; positive metres: supply a known interval",
    )
    parser.add_argument("--height-recovery", choices=["auto", "none"], default="auto")
    parser.add_argument(
        "--height-repair-max-length-m",
        type=float,
        default=10.0,
        help="Maximum total fragment length for the small-piece height repair fallback",
    )
    parser.add_argument(
        "--height-repair-max-fraction",
        type=float,
        default=0.10,
        help="Maximum fraction of a contour loop changed by the small-piece fallback",
    )
    parser.add_argument(
        "--http-client",
        choices=["default", "curl", "httplib"],
        default="default",
        help="DuckDB HTTP client passed to iomaps",
    )
    args = parser.parse_args()
    try:
        run_sheet(
            args.sheet,
            ensure_index(args.index),
            args.output_dir,
            args.cache_dir,
            args.iomaps,
            args.reuse_extracts,
            not args.no_water,
            args.resolution,
            args.spot_radius,
            args.triangle_side,
            args.http_client,
            edge_sliver_m=args.edge_sliver_m,
            form_line_interval=args.form_line_interval,
            contour_contact_radius_m=args.contour_contact_radius_m,
            height_recovery=args.height_recovery,
            height_repair_max_length_m=args.height_repair_max_length_m,
            height_repair_max_fraction=args.height_repair_max_fraction,
            interpolation=args.interpolation,
            local_parquets=args.parquet,
            validate_points=not args.skip_holdout,
        )
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
