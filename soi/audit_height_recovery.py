# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "rasterio>=1.4", "scipy>=1.14", "shapely>=2.1"]
# ///
"""Audit height recovery using a cached classified run; no raster is required.

Writes proposed cycles, a derived contour copy and before/after conflict counts.
Does not change input contours. The DEM CLI uses the same recovery function.
"""

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from pyproj import Transformer
from shapely import segmentize
from shapely.geometry import shape

import contours_to_dem as dem
from contour_contacts import relax_contour_contacts
from contour_height_recovery import recover_contour_heights


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("classified_dir", type=Path)
    parser.add_argument("--sheet", required=True)
    parser.add_argument("--index", type=Path, default=Path("index_50k.geojson"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--resolution", type=float, default=30)
    parser.add_argument("--max-fragment-length-m", type=float, default=10.0)
    parser.add_argument("--max-fragment-fraction", type=float, default=0.10)
    args = parser.parse_args()
    if not math.isfinite(args.resolution) or args.resolution <= 0:
        parser.error("--resolution must be finite and positive")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("Output directory must be empty")
    sheet = dem.get_sheet(args.index, args.sheet)
    epsg = 32600 + min(60, int((sheet.centroid.x + 180) // 6) + 1)
    if sheet.centroid.y < 0:
        epsg += 100
    project = Transformer.from_crs(4326, epsg, always_xy=True).transform
    unproject = Transformer.from_crs(epsg, 4326, always_xy=True).transform
    projected_sheet = dem.transform(project, segmentize(sheet, 0.002))
    classification = json.loads(
        (args.classified_dir / "classification.json").read_text()
    )
    interval = classification.get("sheet_interval") or classification.get(
        "form_line_detection", {}
    ).get("sheet_interval", {})

    def load(name):
        return [
            (dem.transform(project, shape(f["geometry"])), f["properties"])
            for f in dem.read_features(args.classified_dir / name)
        ]

    contours, spots = load("contours.geojsonl"), load("spots.geojsonl")
    derived, audit = recover_contour_heights(
        contours, projected_sheet, interval, boundary_margin_m=args.resolution,
        max_fragment_length_m=args.max_fragment_length_m,
        max_fragment_fraction=args.max_fragment_fraction,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    dem.write_vectors(
        args.output_dir / "recovered_contours.geojsonl", derived, unproject
    )
    dem.write_vectors(args.output_dir / "recovery_cycles.geojsonl", audit, unproject)
    results = {}
    for name, lines in [("before", contours), ("after", derived)]:
        kept, contacts, _ = relax_contour_contacts(
            lines, args.resolution / 2, projected_sheet
        )
        xy, heights = dem.constraint_samples(kept, args.resolution / 2)
        if spots:
            xy = np.concatenate([xy, [[g.x, g.y] for g, _ in spots]])
            heights = np.concatenate([heights, [p["VALUE"] for _, p in spots]])
        conflict_path = args.output_dir / f"{name}_conflicts.geojsonl"
        try:
            dem.unique_samples(
                xy, heights, conflict_path=conflict_path, unproject=unproject
            )
            count = 0
        except ValueError as exc:
            if "conflicting elevations" not in str(exc):
                raise
            count = len(conflict_path.read_text().splitlines())
        results[name] = {
            "height_conflicts": count,
            "cartographic_contacts": len(contacts),
        }
    changes = [p["height_recovery"] for _, p in derived if "height_recovery" in p]
    report = {
        "sheet": args.sheet,
        "scope": "Vector-only recovery on a derived copy, using the same function as the DEM pipeline.",
        "fragment_limits": {"maximum_length_m": args.max_fragment_length_m, "maximum_fraction": args.max_fragment_fraction},
        "algorithm_sha256": hashlib.sha256(Path(__file__).with_name("contour_height_recovery.py").read_bytes()).hexdigest(),
        "regular_contour_interval_m": interval.get("interval_m"),
        "decisions": dict(Counter(p["decision"] for _, p in audit)),
        "changed_source_parts": len(changes),
        "checks": results,
        "point_heights_used_unchanged": len(spots),
        "source_files_changed": False,
        "results": [
            {"coordinates": list(dem.transform(unproject, g).coords), **p}
            for g, p in audit
        ],
        "inputs": [
            {
                "path": str(args.classified_dir / name),
                "sha256": hashlib.sha256(
                    (args.classified_dir / name).read_bytes()
                ).hexdigest(),
            }
            for name in ["contours.geojsonl", "spots.geojsonl", "classification.json"]
        ],
        "limitations": "Requires an isolated closed loop, trustworthy enclosing regular contours, and one inferred sheet interval. Open contours, missing nesting references and inconsistent neighbours remain unresolved. Consistent errors shared by the reference contours cannot be ruled out by this test.",
    }
    (args.output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: v
                for k, v in report.items()
                if k not in {"results", "inputs", "limitations"}
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
