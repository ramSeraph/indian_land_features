# SPDX-License-Identifier: GPL-3.0-or-later
"""Canonical local edge/corner jobs for the SOI production workflow.

Edge jobs always use immutable first-pass DEMs. Corner jobs use those DEMs plus
the same, canonically ordered incident edge patches. Both owners therefore get
the identical shared patch regardless of sheet processing order.
"""

import hashlib
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import rasterio
from affine import Affine
from pyproj import Transformer
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from shapely import STRtree, segmentize
from shapely.geometry import Point, mapping, shape
from shapely.ops import unary_union

from contours_to_dem import transform, write_raster
from seam_interpolate import reconcile_patch


def checksum(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


class SheetIndex:
    def __init__(self, index, selected):
        features = json.loads(Path(index).read_text())["features"]
        self.polygons = {f["properties"]["id"]: shape(f["geometry"]) for f in features}
        selected = sorted(selected)
        if any(s not in self.polygons for s in selected):
            raise ValueError("Selection contains an unknown sheet")
        self.jobs = {}
        self.by_sheet = {s: [] for s in selected}
        polygons = [self.polygons[s] for s in selected]
        tree = STRtree(polygons)
        for i, polygon in enumerate(polygons):
            for j in tree.query(polygon, predicate="intersects"):
                if j <= i:
                    continue
                shared = polygon.boundary.intersection(polygons[j].boundary)
                if shared.length > 0:
                    self._add("edge", [selected[i], selected[j]], shared)
        vertices = {}
        for sheet in selected:
            polygon = self.polygons[sheet]
            for xy in {tuple(c[:2]) for c in polygon.exterior.coords}:
                vertices.setdefault(xy, []).append(sheet)
        for xy, members in sorted(vertices.items()):
            if len(members) >= 3:
                self._add("corner", members, Point(xy))

    def _add(self, kind, members, geometry):
        members = sorted(members)
        identifier = kind + "-" + digest([members, geometry.wkb_hex])[:20]
        self.jobs[identifier] = {
            "id": identifier,
            "kind": kind,
            "members": members,
            "geometry": mapping(geometry),
        }
        for member in members:
            self.by_sheet[member].append(identifier)

    def incident_edges(self, job):
        point = shape(job["geometry"])
        return sorted(
            key
            for key, candidate in self.jobs.items()
            if candidate["kind"] == "edge"
            and set(candidate["members"]).issubset(job["members"])
            and shape(candidate["geometry"]).intersects(point)
        )

    def neighbourhood(self, sheet):
        return sorted(
            {s for key in self.by_sheet[sheet] for s in self.jobs[key]["members"]}
            | {sheet}
        )


def grid_slices(source_shape, source_transform, target_shape, target_transform):
    if (
        source_transform.a != target_transform.a
        or source_transform.e != target_transform.e
    ):
        raise ValueError("Different raster resolutions")
    if any(
        (source_transform.b, source_transform.d, target_transform.b, target_transform.d)
    ):
        raise ValueError("Rotated grids are not supported")
    col, row = ~target_transform * (source_transform.c, source_transform.f)
    if abs(col - round(col)) > 1e-6 or abs(row - round(row)) > 1e-6:
        raise ValueError("Raster origins do not share the same pixel lattice")
    col, row = round(col), round(row)
    r0, r1 = max(0, row), min(target_shape[0], row + source_shape[0])
    c0, c1 = max(0, col), min(target_shape[1], col + source_shape[1])
    if r0 >= r1 or c0 >= c1:
        return None
    return np.s_[r0 - row : r1 - row, c0 - col : c1 - col], np.s_[r0:r1, c0:c1]


def apply_saved_patch(data, affine, crs, inside, patch, changed=None):
    report = json.loads((Path(patch) / "report.json").read_text())
    if crs.to_string() != report["crs"]:
        raise ValueError("Patch and output CRS differ")
    if checksum(Path(patch) / "patch.npz") != report["patch_sha256"]:
        raise ValueError("Saved overlap patch changed")
    with np.load(Path(patch) / "patch.npz") as saved:
        values, mask = saved["values"], saved["mask"]
    placement = grid_slices(
        values.shape, Affine(*report["transform"]), data.shape, affine
    )
    if placement is None:
        return
    source, target = placement
    use = mask[source] & inside[target]
    data[target][use] = values[source][use]
    if changed is not None:
        changed[target][use] = True


def create_patch(
    job,
    folders,
    index,
    output,
    *,
    prior_patches=(),
    width_m=2000,
    context_m=3600,
    fingerprint="",
):
    """Solve one canonical edge or corner, storing the complete numerical inputs."""
    output = Path(output)
    if output.exists():
        raise ValueError("Patch destination already exists")
    members = sorted(folders)
    if (job["kind"] == "edge" and len(members) != 2) or len(members) < 2:
        raise ValueError("A shared edge needs both source sheets")
    metadata = {}
    for sheet in members:
        with rasterio.open(folders[sheet] / "dem.tif") as src:
            metadata[sheet] = (src.shape, src.transform, src.crs)
    _, _, crs = metadata[members[0]]
    resolution = metadata[members[0]][1].a
    if any(m[2] != crs for m in metadata.values()):
        raise ValueError("Cross-UTM overlap is unsupported; retained for a later pass")
    if width_m <= 0 or context_m < 120 * resolution:
        raise ValueError(
            "Patch needs a positive width and a full directional-search context"
        )
    project = Transformer.from_crs(4326, crs, always_xy=True).transform
    domain = transform(
        project, segmentize(unary_union([index.polygons[s] for s in members]), 0.002)
    )
    join = transform(project, segmentize(shape(job["geometry"]), 0.002))
    if job["kind"] == "edge":
        region = join.buffer(width_m)
    else:
        # Cover every intersection of incident strips, including UTM grid
        # convergence. An axis-aligned square can miss their rotated corners.
        bands = [
            transform(
                project, segmentize(shape(index.jobs[key]["geometry"]), 0.002)
            ).buffer(width_m)
            for key in index.incident_edges(job)
            if set(index.jobs[key]["members"]).issubset(members)
        ]
        region = unary_union([a.intersection(b) for a, b in combinations(bands, 2)])
        if region.is_empty:
            raise ValueError("No incident edge overlap is available at this corner")
    region = region.intersection(domain)
    context = region.buffer(context_m).intersection(domain)
    west, south, east, north = context.bounds
    west, south = (
        math.floor(west / resolution) * resolution,
        math.floor(south / resolution) * resolution,
    )
    east, north = (
        math.ceil(east / resolution) * resolution,
        math.ceil(north / resolution) * resolution,
    )
    affine = from_origin(west, north, resolution, resolution)
    grid_shape = (
        round((north - south) / resolution),
        round((east - west) / resolution),
    )
    core = geometry_mask([mapping(region)], grid_shape, affine, invert=True)
    known, prior = np.full(grid_shape, np.nan), np.full(grid_shape, np.nan)
    unknown = np.zeros(grid_shape, dtype=np.int32)
    ownership = np.zeros(grid_shape, dtype=np.int16)
    unknown_mapping, inputs = [], {}
    next_group = 0
    for owner, sheet in enumerate(members, 1):
        folder = folders[sheet]
        old_shape, old_transform, _ = metadata[sheet]
        placement = grid_slices(old_shape, old_transform, grid_shape, affine)
        if placement is None:
            continue
        source, target = placement
        inside = geometry_mask(
            [mapping(transform(project, segmentize(index.polygons[sheet], 0.002)))],
            old_shape,
            old_transform,
            invert=True,
        )
        with rasterio.open(folder / "dem.tif") as src:
            data = src.read(1, masked=True).filled(np.nan).astype(float)
        for patch in prior_patches:
            apply_saved_patch(data, old_transform, crs, inside, patch)
        with np.load(folder / "interpolation_inputs.npz") as saved:
            k = saved["known_heights"].copy()
            u = saved["unknown_ids"].copy()
        k[~inside], u[~inside], data[~inside] = np.nan, 0, np.nan
        full_core = np.zeros(old_shape, dtype=bool)
        full_core[source] = core[target]
        for group in np.unique(u[full_core & (u > 0)]):
            if np.any((u == group) & ~full_core):
                raise ValueError(
                    f"{sheet}: unknown contour {group} crosses the fixed patch boundary"
                )
        # Do not guess which unknown fragments from separate sheets correspond.
        if np.any(u > 0):
            shared = unary_union(
                [
                    index.polygons[sheet].boundary.intersection(
                        index.polygons[s].boundary
                    )
                    for s in members
                    if s != sheet
                ]
            )
            shared = transform(project, segmentize(shared, 0.002))
            near = geometry_mask(
                [mapping(shared.buffer(resolution))],
                old_shape,
                old_transform,
                invert=True,
            )
            if np.any((u > 0) & near):
                raise ValueError(
                    f"{sheet}: an unknown contour touches a shared boundary"
                )
        selected = inside[source]
        if np.any(selected & (ownership[target] != 0)):
            raise ValueError("Sheet interiors overlap on the common grid")
        ownership[target][selected] = owner
        known[target][selected] = k[source][selected]
        prior[target][selected] = data[source][selected]
        for group in np.unique(u[source][u[source] > 0]):
            next_group += 1
            unknown[target][u[source] == group] = next_group
            unknown_mapping.append(
                {"joint_id": next_group, "sheet": sheet, "original_id": int(group)}
            )
        inputs[sheet] = {
            name: checksum(folder / name)
            for name in ("dem.tif", "interpolation_inputs.npz")
        }
    core &= ownership > 0
    patched, numerical = reconcile_patch(known, unknown, prior, core)
    output.mkdir(parents=True)
    np.savez_compressed(
        output / "inputs.npz",
        known=known,
        unknown=unknown,
        prior=prior,
        core=core,
        transform=np.array(tuple(affine)[:6]),
        epsg=crs.to_epsg(),
    )
    rr, cc = np.nonzero(core)
    crop = np.s_[rr.min() : rr.max() + 1, cc.min() : cc.max() + 1]
    patch_transform = affine * Affine.translation(int(cc.min()), int(rr.min()))
    np.savez_compressed(
        output / "patch.npz", values=patched[crop].astype("float32"), mask=core[crop]
    )
    report = {
        **numerical,
        "job": job,
        "members_used": members,
        "run_fingerprint": fingerprint,
        "width_per_side_m": width_m,
        "context_m": context_m,
        "crs": crs.to_string(),
        "transform": list(tuple(patch_transform)[:6]),
        "source_inputs": inputs,
        "unknown_id_mapping": unknown_mapping,
        "prior_edge_patches": [
            {
                "job": json.loads((p / "report.json").read_text())["job"]["id"],
                "report_sha256": checksum(p / "report.json"),
            }
            for p in prior_patches
        ],
        "inputs_sha256": checksum(output / "inputs.npz"),
        "patch_sha256": checksum(output / "patch.npz"),
    }
    (output / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


def correct_sheet(sheet, directory, index, patches, outcomes, *, width_m=2000):
    """Apply shared artifacts to a copied first-pass sheet and validate invariants."""
    directory = Path(directory)
    path = directory / "dem.tif"
    with rasterio.open(path) as src:
        data = src.read(1, masked=True).filled(np.nan)
        affine, crs = src.transform, src.crs
    original = data.copy()
    project = Transformer.from_crs(4326, crs, always_xy=True).transform
    inside = geometry_mask(
        [mapping(transform(project, segmentize(index.polygons[sheet], 0.002)))],
        data.shape,
        affine,
        invert=True,
    )
    changed = np.zeros(data.shape, dtype=bool)
    ordered = sorted(
        patches, key=lambda key: (index.jobs[key]["kind"] == "corner", key)
    )
    for key in ordered:
        apply_saved_patch(data, affine, crs, inside, patches[key], changed)
    with np.load(directory / "interpolation_inputs.npz") as inputs:
        known = inputs["known_heights"]
    measured = np.isfinite(known) & inside
    if not np.array_equal(data[measured], known[measured].astype("float32")):
        raise ValueError("An overlap patch changed a measured height")
    if not np.array_equal(data[~changed], original[~changed], equal_nan=True):
        raise ValueError("An overlap patch changed a cell outside its mask")
    if np.any(np.isfinite(original) & ~np.isfinite(data)):
        raise ValueError("An overlap patch removed previously supported data")
    if np.any(np.isfinite(data) & ~inside):
        raise ValueError("An overlap patch escaped the sheet mask")
    unresolved = {
        k: v
        for k, v in outcomes.items()
        if v["status"] != "applied" or v.get("missing_members")
    }
    status = (
        "partial"
        if unresolved
        else (
            "complete_for_available_neighbours" if ordered else "no_selected_neighbours"
        )
    )
    summary = {
        "status": status,
        "seam_reconciled": not unresolved and bool(ordered),
        "width_per_side_m": width_m,
        "jobs": outcomes,
        "applied_order": ordered,
        "fixed_constraint_pixels_checked": int(measured.sum()),
        "outside_patches_unchanged": True,
        "missing_coverage_is_not_filled": True,
    }
    # Keep the immutable first pass in the audit package, with its original report.
    path.replace(directory / "dem_initial.tif")
    (directory / "report.json").replace(directory / "report_initial.json")
    write_raster(
        path,
        data,
        affine,
        crs.to_epsg(),
        "SOI elevation with local shared-edge correction",
    )
    with rasterio.open(path, "r+") as saved:
        saved.update_tags(
            SOI_SEAM_STATUS=status,
            SOI_SEAM_WIDTH_M=str(width_m),
            SOI_WATER_CONDITIONING="none",
        )
    report = json.loads((directory / "report_initial.json").read_text())
    report.update(
        seam_reconciled=summary["seam_reconciled"],
        overlap=summary,
        valid_cells=int(np.isfinite(data).sum()),
        coverage_fraction=float(np.isfinite(data).sum() / inside.sum()),
    )
    (directory / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    (directory / "overlap.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    return summary
