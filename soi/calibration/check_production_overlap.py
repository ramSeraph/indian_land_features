# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "rasterio>=1.4", "scipy>=1.14", "pyproj>=3.6", "shapely>=2.1"]
# ///
"""Check packaged overlap outputs against GLO-30 and replay archived patches."""

import argparse
import io
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.features import geometry_mask
from rasterio.io import MemoryFile
from shapely import segmentize
from shapely.geometry import mapping
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from contours_to_dem import transform
from seam_interpolate import reconcile_patch
from sheet_overlaps import SheetIndex, checksum


def metrics(values):
    values = np.asarray(values, dtype=float)
    return {
        "count": len(values),
        "mean_m": float(values.mean()),
        "rmse_m": float(np.sqrt(np.mean(values**2))),
    }


def check(root, reference):
    root, reference = Path(root), Path(reference)
    state = json.loads((root / "state.json").read_text())
    selected = json.loads((root / "run.json").read_text())["configuration"]["sheets"]
    index = SheetIndex(ROOT / "index_50k.geojson", selected)
    results, replayed = {}, {}
    before_errors, after_errors = [], []
    for sheet, record in state["sheets"].items():
        if record["outcome"]["status"] != "built":
            continue
        folder = root / "packages" / sheet
        base = f"SOI_DEM_30m_{sheet}"
        with zipfile.ZipFile(folder / f"{base}.audit.zip") as archive:
            assert archive.testzip() is None
            overlap = json.loads(archive.read("dem/overlap.json"))
            with (
                MemoryFile(archive.read("dem/dem_initial.tif")) as memory,
                memory.open() as old,
            ):
                before = old.read(1, masked=True).filled(np.nan)
                old_transform, old_crs = old.transform, old.crs
            with rasterio.open(folder / f"{base}.tif") as final:
                after = final.read(1, masked=True).filled(np.nan)
                assert final.transform == old_transform and final.crs == old_crs
                assert final.tags(ns="IMAGE_STRUCTURE")["LAYOUT"] == "COG"
                assert final.tags()["SOI_SEAM_STATUS"] == overlap["status"]
                affine, crs = final.transform, final.crs
            with rasterio.open(reference / sheet / "glo30_aligned.tif") as glo:
                assert (
                    glo.transform == affine
                    and glo.crs == crs
                    and glo.shape == after.shape
                )
                heights = glo.read(1, masked=True).filled(np.nan)
            edges = [
                index.jobs[key]
                for key in index.by_sheet[sheet]
                if index.jobs[key]["kind"] == "edge"
            ]
            joins = unary_union(
                [
                    index.polygons[a].boundary.intersection(index.polygons[b].boundary)
                    for edge in edges
                    for a, b in [edge["members"]]
                ]
            )
            project = Transformer.from_crs(4326, crs, always_xy=True).transform
            near = geometry_mask(
                [mapping(transform(project, segmentize(joins, 0.002)).buffer(300))],
                after.shape,
                affine,
                invert=True,
            )
            common = (
                near & np.isfinite(before) & np.isfinite(after) & np.isfinite(heights)
            )
            a, b = (before - heights)[common], (after - heights)[common]
            before_errors.extend(a)
            after_errors.extend(b)
            results[sheet] = {
                "before": metrics(a),
                "after": metrics(b),
                "overlap": overlap["status"],
                "height_points_checked": record["outcome"]["checks"][
                    "height_points_checked"
                ],
                "cog_sha256": checksum(folder / f"{base}.tif"),
                "archive_sha256": checksum(folder / f"{base}.audit.zip"),
            }
            for name in archive.namelist():
                if not name.startswith("overlaps/") or not name.endswith("/inputs.npz"):
                    continue
                prefix = name.removesuffix("inputs.npz")
                report = json.loads(archive.read(prefix + "report.json"))
                key = report["job"]["id"]
                if key in replayed:
                    assert replayed[key]["patch_sha256"] == report["patch_sha256"]
                    continue
                with np.load(io.BytesIO(archive.read(name))) as saved:
                    known, unknown, prior, core = [
                        saved[n] for n in ("known", "unknown", "prior", "core")
                    ]
                restored, _ = reconcile_patch(known, unknown, prior, core)
                rr, cc = np.nonzero(core)
                cropped = restored[
                    rr.min() : rr.max() + 1, cc.min() : cc.max() + 1
                ].astype("float32")
                with np.load(io.BytesIO(archive.read(prefix + "patch.npz"))) as patch:
                    np.testing.assert_array_equal(cropped, patch["values"])
                replayed[key] = {
                    "status": "pixel_identical",
                    "patch_sha256": report["patch_sha256"],
                }
    result = {
        "sheets": results,
        "within_300m_of_shared_edges": {
            "before": metrics(before_errors),
            "after": metrics(after_errors),
        },
        "archived_patch_replay": replayed,
        "script_sha256": checksum(__file__),
        "run_manifest_sha256": checksum(root / "run.json"),
    }
    (root / "overlap_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument(
        "--reference-dir", type=Path, default=Path("data/comparisons/production-pilot")
    )
    args = parser.parse_args()
    check(args.run, args.reference_dir)
