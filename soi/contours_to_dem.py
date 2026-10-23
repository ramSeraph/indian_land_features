#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy>=2", "pyproj>=3.6", "rasterio>=1.4", "scipy>=1.14", "shapely>=2.1"]
# ///
"""Build a single-sheet DEM from an iomaps GeoJSONL extract.

The default is a Python adaptation of Imagico's contour interpolation stages.
The previous linear TIN remains available explicitly for baseline comparisons.
Input coordinates must be WGS84 longitude/latitude; VALUE is elevation in metres.
"""

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree
from shapely import STRtree, line_merge, normalize, segmentize
from shapely import transform as transform_coordinates
from shapely.geometry import LineString, MultiLineString, Point, Polygon, mapping, shape
from shapely.ops import unary_union

from contour_form_lines import (
    contour_chains,
    identify_form_line_candidates,
    infer_sheet_interval,
    parse_form_line_interval,
)
from contour_contacts import relax_contour_contacts
from contour_height_recovery import recover_contour_heights, validate_recovery_limits
from imagico_dem import build_grid as build_imagico_grid, sample_grid

NODATA = -9999.0


def transform(function, geometry):
    return transform_coordinates(geometry, function, interleaved=False)


def read_features(path):
    """Read GeoJSONSeq, including optional RFC 8142 record separators."""
    with Path(path).open() as stream:
        for number, line in enumerate(stream, 1):
            line = line.strip().lstrip("\x1e")
            if line:
                feature = json.loads(line)
                if feature.get("type") != "Feature":
                    raise ValueError(f"{path}:{number}: expected a GeoJSON Feature")
                yield feature


def get_sheet(index, sheet_id):
    features = json.loads(Path(index).read_text())["features"]
    matches = [f for f in features if f["properties"].get("id") == sheet_id]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one sheet {sheet_id!r}; found {len(matches)}"
        )
    return shape(matches[0]["geometry"])


def parts(geometry, kind):
    if geometry.geom_type == kind:
        yield geometry
    elif (
        geometry.geom_type.startswith("Multi")
        or geometry.geom_type == "GeometryCollection"
    ):
        for part in geometry.geoms:
            yield from parts(part, kind)


def roundness(line):
    """Isoperimetric quotient of a simple closed line; a circle approaches 1."""
    if line.geom_type != "LineString" or not line.is_ring or line.length == 0:
        return 0.0
    return float(4 * math.pi * Polygon(line).area / line.length**2)


def geometry_hash(geometry):
    return hashlib.sha256(normalize(geometry).wkb).hexdigest()


def spot_ring(line, max_radius_m=13.0, tolerance=0.03, min_roundness=0.98):
    """Identify small closed circular symbols in projected metre coordinates.

    A geometry-only heuristic: return the fitted centre and radius for audit.
    Real small circular contours can also pass, so retain original source data.
    """
    if line.geom_type != "LineString" or not line.is_ring or len(line.coords) < 9:
        return None
    # Points sampled along an arc can fit a circle perfectly even when the
    # closing chord makes the polygon decidedly noncircular.
    if roundness(line) < min_roundness:
        return None
    xy = np.asarray(line.coords)[:-1, :2]
    origin = xy.mean(axis=0)
    local = xy - origin
    system = np.column_stack((2 * local, np.ones(len(local))))
    solution, _, rank, _ = np.linalg.lstsq(system, np.sum(local**2, axis=1), rcond=None)
    if rank < 3:
        return None
    centre = origin + solution[:2]
    distances = np.linalg.norm(xy - centre, axis=1)
    radius = float(distances.mean())
    if 0 < radius <= max_radius_m and distances.std() / radius <= tolerance:
        return Point(centre), radius
    return None


def spot_triangle(line, nominal_side_m=70.0):
    """Find closed triangular height symbols at the calibrated map-symbol scale.

    Sides must be within 20% of the nominal length and the shortest at least
    80% of the longest. Allow 1 m outline error at the default 70 m scale, so
    collinear or slightly noisy extra vertices do not hide the three corners.
    A real contour with the same geometry can still pass this heuristic.
    """
    if (
        nominal_side_m <= 0
        or line.geom_type != "LineString"
        or not line.is_ring
        or not line.is_valid
        or line.length == 0
    ):
        return None
    edge_tolerance = nominal_side_m / 70.0
    max_side = nominal_side_m * 1.2
    minx, miny, maxx, maxy = line.bounds
    if max(maxx - minx, maxy - miny) > max_side + 2 * edge_tolerance:
        return None
    polygon = Polygon(line)
    # The hull removes collinear samples without relying on the ring's starting
    # vertex. The outline/area checks below reject concave or curved impostors.
    # Topology-preserving simplification can retain a redundant fourth corner.
    # A convex hull needs no such protection; validate the resulting type below.
    triangle = polygon.convex_hull.simplify(edge_tolerance, preserve_topology=False)
    if triangle.geom_type != "Polygon" or len(triangle.exterior.coords) != 4:
        return None
    sides = np.linalg.norm(
        np.diff(np.asarray(triangle.exterior.coords)[:, :2], axis=0), axis=1
    )
    shortest, longest = float(sides.min()), float(sides.max())
    if (
        shortest < nominal_side_m * 0.8
        or longest > max_side
        or shortest / longest < 0.8
    ):
        return None
    # Densification checks the closing segments too, not just original vertices.
    edge_error = float(
        segmentize(line, edge_tolerance).hausdorff_distance(
            segmentize(triangle.exterior, edge_tolerance)
        )
    )
    area_error = abs(polygon.area / triangle.area - 1)
    if edge_error > edge_tolerance or area_error > 0.03:
        return None
    return polygon.centroid, {
        "symbol_side_min_m": shortest,
        "symbol_side_max_m": longest,
        "symbol_edge_error_m": edge_error,
        "symbol_area_error_fraction": float(area_error),
    }


def reconstruct_spot_symbols(
    features, project, max_radius_m=13.0, triangle_side_m=70.0
):
    """Recover split symbols using exact endpoints and identical elevations.

    Only complete rings passing the existing symbol tests are reconstructed.
    No snapping, missing-edge inference, or elevation correction is performed.
    Every contributing source component retains its attributes and geometry hash.
    """
    groups, seen = defaultdict(list), set()
    max_length = max(2 * math.pi * max_radius_m, 3 * 1.2 * triangle_side_m)
    for feature in features:
        properties = feature.get("properties") or {}
        try:
            elevation = float(properties["VALUE"])
        except (KeyError, ValueError, TypeError):
            continue  # The normal validation pass reports malformed records.
        geometry = shape(feature["geometry"])
        for part_index, component in enumerate(parts(geometry, "LineString")):
            if component.is_closed or not component.is_valid:
                continue
            key = normalize(component).wkb, elevation
            if key in seen:
                continue
            seen.add(key)
            if transform(project, component).length > max_length:
                continue
            groups[elevation].append(
                (
                    component,
                    key,
                    {
                        "source_id": properties.get("id", feature.get("id")),
                        "source_part_index": part_index,
                        "source_geometry_sha256": geometry_hash(component),
                        "properties": properties,
                    },
                )
            )
    recovered = {}
    for records in groups.values():
        if len(records) < 2:
            continue
        geometries = [r[0] for r in records]
        tree = STRtree(geometries)
        for ring in parts(line_merge(MultiLineString(geometries)), "LineString"):
            if not ring.is_ring:
                continue
            projected = transform(project, ring)
            if (
                spot_ring(projected, max_radius_m) is None
                and spot_triangle(projected, triangle_side_m) is None
            ):
                continue
            indices = sorted(int(i) for i in tree.query(ring, predicate="covers"))
            if len(indices) < 2:
                continue
            provenance = [records[i][2] for i in indices]
            for i in indices:
                recovered[records[i][1]] = ring, provenance
    return recovered


def prepare_constraints(
    paths,
    sheet,
    project,
    max_radius_m=13.0,
    reviewed_spots=(),
    *,
    triangle_side_m=70.0,
    sheet_id=None,
    edge_sliver_m=0.0,
    form_line_interval="auto",
    excluded_contours=None,
    form_line_audit=None,
    form_line_interval_audit=None,
):
    """Classify before clipping; optionally collect complete excluded line parts.

    Edge slivers are mostly outside this sheet and have no retained geometry
    beyond the given inward boundary band. Removing them is a local extraction
    policy, not a reconciliation of heights with the neighbouring sheet.
    """
    if not math.isfinite(edge_sliver_m) or edge_sliver_m < 0:
        raise ValueError("Edge sliver tolerance must be nonnegative and finite")
    contours, spots, raw_contours = [], [], []
    counts, elevations = Counter(), Counter()
    seen = set()
    reviews = {review["geometry_sha256"]: review for review in reviewed_spots}
    if len(reviews) != len(reviewed_spots):
        raise ValueError("Duplicate geometry in reviewed spot corrections")
    used_reviews = set()
    projected_sheet = transform(project, segmentize(sheet, 0.002))
    reconstructed = reconstruct_spot_symbols(
        (feature for path in paths for feature in read_features(path)),
        project,
        max_radius_m,
        triangle_side_m,
    )
    used_reconstructions = set()
    for path in paths:
        for feature in read_features(path):
            counts["input_features"] += 1
            properties = feature.get("properties") or {}
            try:
                elevation = float(properties["VALUE"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid or missing VALUE in {path}: {properties}"
                ) from exc
            if not math.isfinite(elevation) or elevation == NODATA:
                raise ValueError(f"Invalid elevation: {elevation}")
            geometry = shape(feature["geometry"])
            if geometry.geom_type not in ("LineString", "MultiLineString", "Point"):
                raise ValueError(f"Unexpected geometry: {geometry.geom_type}")
            components = (
                geometry.geoms
                if geometry.geom_type == "MultiLineString"
                else [geometry]
            )
            original_properties = properties
            for part_index, component in enumerate(components):
                properties = original_properties
                key = (normalize(component).wkb, elevation)
                if key in seen:
                    counts["duplicate_parts"] += 1
                    continue
                seen.add(key)
                source_fragments = None
                if key in reconstructed:
                    component, source_fragments = reconstructed[key]
                    joined_key = geometry_hash(component), elevation
                    if joined_key in used_reconstructions:
                        continue
                    used_reconstructions.add(joined_key)
                    properties = source_fragments[0]["properties"]
                    part_index = source_fragments[0]["source_part_index"]
                projected = transform(project, component)
                if projected.geom_type == "LineString" and projected.length == 0:
                    counts["zero_length_line_parts"] += 1
                    continue
                ring = spot_ring(projected, max_radius_m)
                triangle = (
                    spot_triangle(projected, triangle_side_m) if ring is None else None
                )
                metadata = dict(properties, VALUE=elevation)
                review_key = geometry_hash(component)
                review = reviews.get(review_key)
                if review:
                    if (
                        properties.get("id") != review["source_id"]
                        or abs(elevation - float(review["value_m"])) > 1e-9
                        or projected.geom_type != "LineString"
                        or not projected.is_ring
                    ):
                        raise ValueError(
                            f"Reviewed spot no longer matches source: {review['source_id']}"
                        )
                    used_reviews.add(review_key)
                    point = Polygon(projected).centroid
                    metadata.update(
                        spot_classification="reviewed_symbol",
                        review_note=review["reason"],
                    )
                elif triangle:
                    point, measurements = triangle
                    metadata.update(
                        spot_classification="triangular_symbol", **measurements
                    )
                elif ring is not None or projected.geom_type == "Point":
                    point = ring[0] if ring else projected
                    metadata["spot_classification"] = (
                        "circular_symbol" if ring else "source_point"
                    )
                else:
                    point = None
                if point is not None:
                    if projected_sheet.covers(point):
                        symbol_type = (
                            "triangle"
                            if triangle
                            else "circle"
                            if ring
                            else "point"
                            if projected.geom_type == "Point"
                            else "unknown"
                        )
                        metadata.update(
                            symbol_type=symbol_type,
                            elevation_m=elevation,
                            source_id=properties.get("id", feature.get("id")),
                            source_part_index=part_index,
                            source_geometry_sha256=review_key,
                            classification_method=(
                                "manual_review"
                                if review
                                else "source_geometry"
                                if symbol_type == "point"
                                else "geometry_heuristic"
                            ),
                        )
                        if sheet_id is not None:
                            metadata["source_sheet"] = sheet_id
                        if source_fragments:
                            metadata["source_fragments"] = source_fragments
                            metadata["classification_method"] = (
                                "geometry_heuristic_joined_fragments"
                            )
                            counts["reconstructed_spot_heights"] += 1
                            counts["reconstructed_symbol_fragments"] += len(
                                source_fragments
                            )
                        if triangle:
                            metadata.update(triangle[1])
                        if ring:
                            metadata["symbol_radius_m"] = ring[1]
                            metadata["symbol_roundness"] = roundness(projected)
                        if review:
                            counts["reviewed_spot_heights"] += 1
                        elif triangle:
                            counts["triangular_spot_heights"] += 1
                        elif ring:
                            counts["circular_spot_heights"] += 1
                        spots.append((point, metadata))
                else:
                    metadata.update(
                        source_id=properties.get("id", feature.get("id")),
                        source_part_index=part_index,
                        source_geometry_sha256=review_key,
                    )
                    if sheet_id is not None:
                        metadata["source_sheet"] = sheet_id
                    raw_contours.append((component, projected, metadata))
    if used_reviews != reviews.keys():
        raise ValueError(
            "Reviewed spot geometry missing from input; check the source version"
        )
    form_lines = identify_form_line_candidates(
        raw_contours,
        form_line_interval,
        audit=form_line_audit,
        sheet=projected_sheet,
        interval_audit=form_line_interval_audit,
    )
    interior = projected_sheet.buffer(-edge_sliver_m)
    if edge_sliver_m and interior.is_empty:
        raise ValueError("Edge sliver tolerance leaves no sheet interior")
    for number, (component, projected, metadata) in enumerate(raw_contours):
        # Clip in the source CRS. Projected straight chords can fall just outside
        # the projected curved boundary and lose lines.
        clipped = component.intersection(sheet)
        retained = transform(project, segmentize(clipped, 0.002))
        if retained.length == 0:
            continue
        reason, evidence = None, {}
        if number in form_lines:
            reason = "form_line_candidate"
            evidence = form_lines[number]
        elif (
            edge_sliver_m > 0
            and retained.length < projected.length / 2
            and not retained.intersects(interior)
        ):
            reason = "outside_edge_sliver"
            evidence = {
                "edge_tolerance_m": edge_sliver_m,
                "source_length_m": projected.length,
                "retained_length_m": retained.length,
                "retained_fraction": retained.length / projected.length,
                "seam_reconciled": False,
            }
        if reason:
            counts[f"excluded_{reason}_parts"] += 1
            if excluded_contours is not None:
                excluded_contours.append(
                    (
                        projected,
                        dict(
                            metadata,
                            exclusion_reason=reason,
                            exclusion_evidence=evidence,
                            used_in_dem=False,
                        ),
                    )
                )
            continue
        for line in parts(clipped, "LineString"):
            if line.length > 0:
                contours.append((transform(project, segmentize(line, 0.002)), metadata))
                elevations[metadata["VALUE"]] += 1
    counts.update(contour_parts=len(contours), spot_heights=len(spots))
    return (
        contours,
        spots,
        projected_sheet,
        dict(counts),
        dict(sorted(elevations.items())),
    )


def constraint_samples(contours, spacing):
    samples, heights = [], []
    for line, properties in contours:
        xy = np.asarray(segmentize(line, spacing).coords)[:, :2]
        samples.append(xy)
        heights.append(np.full(len(xy), properties["VALUE"]))
    if not samples:
        raise ValueError("No contours intersect the requested sheet")
    return np.concatenate(samples), np.concatenate(heights)


def unique_samples(xy, heights, tolerance=0.01, *, conflict_path=None, unproject=None):
    """Coalesce samples within millimetre rounding; reject conflicting heights."""
    _, first, inverse = np.unique(
        np.round(xy, 3), axis=0, return_index=True, return_inverse=True
    )
    low = np.full(len(first), np.inf)
    high = np.full(len(first), -np.inf)
    np.minimum.at(low, inverse, heights)
    np.maximum.at(high, inverse, heights)
    conflicts = int(np.count_nonzero(high - low > tolerance))
    if conflicts:
        if conflict_path is not None:
            records = [
                (
                    Point(xy[first[group]]),
                    {
                        "elevations_m": np.unique(heights[inverse == group]).tolist(),
                        "range_m": float(high[group] - low[group]),
                    },
                )
                for group in np.flatnonzero(high - low > tolerance)
            ]
            write_vectors(conflict_path, records, unproject)
        raise ValueError(
            f"{conflicts} coincident sample locations have conflicting elevations"
        )
    return xy[first], (low + high) / 2


def error_stats(errors):
    errors = np.asarray(errors)
    valid = errors[np.isfinite(errors)]
    result = {"count": int(valid.size), "unavailable": int(errors.size - valid.size)}
    if valid.size:
        result.update(
            bias_m=float(valid.mean()),
            mae_m=float(np.abs(valid).mean()),
            rmse_m=float(np.sqrt(np.mean(valid**2))),
            max_abs_m=float(np.abs(valid).max()),
        )
    return result


def flatten_water(
    dem, water_shapes, model, grid_transform, sheet, resolution, max_shore_range=5.0
):
    """Flatten enclosed polygons to a modelled shoreline median; log every decision.

    Crossing-sheet polygons and shorelines with >5 m variation are left alone.
    No elevations or vertical datum are inferred from CartoDEM attributes.
    """
    result, decisions = dem.copy(), []
    for polygon in parts(unary_union(water_shapes), "Polygon"):
        decision = {"area_m2": float(polygon.area)}
        decisions.append(decision)
        if not sheet.covers(polygon):
            decision["status"] = "skipped_crosses_sheet"
            continue
        mask = geometry_mask([mapping(polygon)], dem.shape, grid_transform, invert=True)
        if not np.any(mask & np.isfinite(dem)):
            decision["status"] = "skipped_no_supported_cells"
            continue
        shoreline = np.concatenate(
            [
                np.asarray(segmentize(LineString(ring.coords), resolution).coords)[
                    :, :2
                ]
                for ring in [polygon.exterior, *polygon.interiors]
            ]
        )
        heights = model(shoreline)
        if not np.all(np.isfinite(heights)):
            decision["status"] = "skipped_unsupported_shoreline"
            continue
        spread = float(np.ptp(heights))
        decision["shore_range_m"] = spread
        if spread > max_shore_range:
            decision["status"] = "skipped_variable_shoreline"
            continue
        level = float(np.median(heights))
        result[mask & np.isfinite(dem)] = level
        decision.update(
            status="flattened",
            elevation_m=level,
            method="interpolated_shoreline_median",
        )
    return result, decisions


def write_raster(path, data, grid_transform, epsg, description, method=None):
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=data.shape[1],
        height=data.shape[0],
        count=1,
        dtype="float32",
        crs=f"EPSG:{epsg}",
        transform=grid_transform,
        nodata=NODATA,
        compress="deflate",
        predictor=3,
        tiled=True,
    ) as dst:
        dst.write(np.where(np.isfinite(data), data, NODATA).astype("float32"), 1)
        dst.set_band_description(1, description)
        dst.set_band_unit(1, "m")
        dst.update_tags(elevation_source="Survey of India VALUE")
        if method:
            dst.update_tags(method=method)


def write_vectors(path, records, unproject):
    with path.open("w") as stream:
        for geometry, properties in records:
            json.dump(
                {
                    "type": "Feature",
                    "properties": properties,
                    "geometry": mapping(transform(unproject, geometry)),
                },
                stream,
            )
            stream.write("\n")


def build_dem(
    inputs,
    index,
    sheet_id,
    output_dir,
    resolution=30.0,
    water=None,
    spot_radius=13.0,
    spot_review=None,
    triangle_side=70.0,
    edge_sliver_m=2.0,
    form_line_interval="auto",
    contour_contact_radius_m=None,
    height_recovery="auto",
    height_repair_max_length_m=10.0,
    height_repair_max_fraction=0.10,
    interpolation="imagico",
    validate_points=True,
):
    output_dir = Path(output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("Output directory is not empty; choose a new run directory")
    if not math.isfinite(resolution) or resolution <= 0:
        raise ValueError("Resolution must be positive and finite")
    if height_recovery not in {"auto", "none"}:
        raise ValueError("Height recovery must be auto or none")
    if interpolation not in {"imagico", "linear"}:
        raise ValueError("Interpolation must be imagico or linear")
    validate_recovery_limits(height_repair_max_length_m, height_repair_max_fraction)
    contact_radius = (
        resolution / 2 if contour_contact_radius_m is None else contour_contact_radius_m
    )
    if not math.isfinite(contact_radius) or contact_radius < 0:
        raise ValueError("Contour contact radius must be nonnegative and finite")
    if not math.isfinite(spot_radius) or spot_radius < 0:
        raise ValueError(
            "Spot radius must be nonnegative and finite; zero disables ring conversion"
        )
    if not math.isfinite(triangle_side) or triangle_side < 0:
        raise ValueError(
            "Triangle side must be nonnegative and finite; zero disables triangle conversion"
        )
    sheet = get_sheet(index, sheet_id)
    zone = min(60, int((sheet.centroid.x + 180) // 6) + 1)
    epsg = (32600 if sheet.centroid.y >= 0 else 32700) + zone
    project = Transformer.from_crs(4326, epsg, always_xy=True).transform
    unproject = Transformer.from_crs(epsg, 4326, always_xy=True).transform
    review_data = json.loads(Path(spot_review).read_text()) if spot_review else {}
    if review_data and review_data.get("sheet") != sheet_id:
        raise ValueError("Spot review is for a different sheet")
    form_line_interval = parse_form_line_interval(form_line_interval)
    excluded_contours, form_line_audit, form_line_interval_audit = [], [], {}
    contours, spots, projected_sheet, counts, levels = prepare_constraints(
        inputs,
        sheet,
        project,
        spot_radius,
        review_data.get("spots", []),
        triangle_side_m=triangle_side,
        sheet_id=sheet_id,
        edge_sliver_m=edge_sliver_m,
        form_line_interval=form_line_interval,
        excluded_contours=excluded_contours,
        form_line_audit=form_line_audit,
        form_line_interval_audit=form_line_interval_audit,
    )
    # Save classifications even when contradictory constraints prevent a DEM.
    # The point dataset is independently useful and must remain reviewable.
    output_dir.mkdir(parents=True, exist_ok=True)
    write_vectors(output_dir / "contours.geojsonl", contours, unproject)
    write_vectors(output_dir / "spots.geojsonl", spots, unproject)
    write_vectors(
        output_dir / "excluded_contours.geojsonl", excluded_contours, unproject
    )
    write_vectors(
        output_dir / "form_line_candidates.geojsonl", form_line_audit, unproject
    )
    form_line_detection = {
        "mode": "auto"
        if form_line_interval == "auto"
        else "disabled"
        if form_line_interval is None
        else "supplied",
        "audit": "form_line_candidates.geojsonl",
        "sheet_interval": form_line_interval_audit,
        "decisions_by_source_part": dict(
            Counter(p["form_line_evidence"]["decision"] for _, p in form_line_audit)
        ),
    }
    triangle_points = [
        (point, properties)
        for point, properties in spots
        if properties["symbol_type"] == "triangle"
    ]
    write_vectors(output_dir / "triangle_points.geojsonl", triangle_points, unproject)
    recovered_contours, height_audit = contours, []
    height_interval = form_line_interval_audit
    if height_recovery == "auto":
        if height_interval.get("status") in {"disabled", "supplied"}:
            height_interval = infer_sheet_interval(
                contour_chains([(g, g, p) for g, p in contours]), {}, projected_sheet
            )
        recovered_contours, height_audit = recover_contour_heights(
            contours,
            projected_sheet,
            height_interval,
            boundary_margin_m=resolution,
            max_fragment_length_m=height_repair_max_length_m,
            max_fragment_fraction=height_repair_max_fraction,
        )
    unknown_contours = []
    if interpolation == "imagico":
        groups = {}
        prepared = []
        for line, props in recovered_contours:
            if "height_recovery" in props:
                selection = props["height_recovery"]
                group = groups.setdefault(selection["recovery_id"], len(groups) + 1)
                props = dict(
                    props,
                    VALUE=None,
                    unknown_height_id=group,
                    original_elevation_m=selection["original_elevation_m"],
                    height_status="unassigned",
                    height_unassignment=dict(selection),
                )
                unknown_contours.append((line, props))
            prepared.append((line, props))
        recovered_contours = prepared
        for _, evidence in height_audit:
            if evidence["decision"] == "supported_recovery":
                evidence["application"] = "unassign_height_for_imagico_interpolation"
    write_vectors(
        output_dir / "unknown_height_contours.geojsonl", unknown_contours, unproject
    )
    write_vectors(
        output_dir / "recovered_contours.geojsonl", recovered_contours, unproject
    )
    write_vectors(output_dir / "height_recovery.geojsonl", height_audit, unproject)
    height_handling = {
        "mode": height_recovery,
        "application": "unassign_then_interpolate"
        if interpolation == "imagico"
        else "direct_vector_repair",
        "interval_m": height_interval.get("interval_m"),
        "interval_status": height_interval.get("status"),
        "supported_cycles": sum(
            p["decision"] == "supported_recovery" for _, p in height_audit
        ),
        "methods": dict(
            Counter(
                p["method"]
                for _, p in height_audit
                if p["decision"] == "supported_recovery"
            )
        ),
        "changed_source_parts": sum(
            a[1]["VALUE"] != b[1]["VALUE"] for a, b in zip(contours, recovered_contours)
        ),
        "fragment_max_length_m": height_repair_max_length_m,
        "fragment_max_fraction": height_repair_max_fraction,
        "boundary_margin_m": resolution,
        "audit": "height_recovery.geojsonl",
        "derived_contours": "recovered_contours.geojsonl",
        "source_contours": "contours.geojsonl",
        "source_files_changed": False,
    }
    interpolation_contours, contacts, relaxed_sections = relax_contour_contacts(
        [(g, p) for g, p in recovered_contours if p["VALUE"] is not None],
        contact_radius,
        projected_sheet,
    )
    write_vectors(output_dir / "contour_contacts.geojsonl", contacts, unproject)
    write_vectors(
        output_dir / "relaxed_contour_sections.geojsonl", relaxed_sections, unproject
    )
    contact_handling = {
        "radius_m": contact_radius,
        "radius_source": "half_grid_spacing"
        if contour_contact_radius_m is None
        else "supplied",
        "contact_locations": len(contacts),
        "geometry_patterns": dict(Counter(p["geometry_pattern"] for _, p in contacts)),
        "relaxed_length_m": sum(line.length for line, _ in relaxed_sections),
        "audit": "contour_contacts.geojsonl",
        "relaxed_sections": "relaxed_contour_sections.geojsonl",
        "source_heights_changed": False,
    }
    (output_dir / "classification.json").write_text(
        json.dumps(
            {
                "sheet": sheet_id,
                "counts": counts,
                "symbol_type_counts": dict(Counter(p["symbol_type"] for _, p in spots)),
                "spot_radius_threshold_m": spot_radius,
                "triangle_nominal_side_m": triangle_side,
                "edge_sliver_tolerance_m": edge_sliver_m,
                "form_line_interval_m": form_line_interval_audit.get("interval_m")
                if form_line_interval == "auto"
                else form_line_interval,
                "form_line_detection": form_line_detection,
                "contour_contact_handling": contact_handling,
                "height_recovery": height_handling,
                "interpolation": interpolation,
                "seam_reconciled": False,
                "source_files": [
                    {
                        "path": str(p),
                        "sha256": hashlib.sha256(Path(p).read_bytes()).hexdigest(),
                    }
                    for p in inputs
                ],
                "note": "Classification outputs; report.json and dem.tif are written only after successful interpolation.",
            },
            indent=2,
        )
        + "\n"
    )
    conflict_options = {
        "conflict_path": output_dir / "conflicting_constraints.geojsonl",
        "unproject": unproject,
    }
    xy, heights = constraint_samples(interpolation_contours, resolution / 2)
    xy, heights = unique_samples(xy, heights, **conflict_options)
    print(
        f"{len(contours)} contour parts, {len(spots)} spot candidates, {len(xy)} contour samples",
        flush=True,
    )
    contour_model = (
        LinearNDInterpolator(xy, heights) if interpolation == "linear" else None
    )
    held_out = {
        "count": 0,
        "unavailable": len(spots),
        "status": "not_run_final_uses_all_points",
    }
    if spots:
        spot_xy = np.array([[g.x, g.y] for g, _ in spots])
        spot_z = np.array([p["VALUE"] for _, p in spots])
        if contour_model is not None and validate_points:
            held_out = error_stats(contour_model(spot_xy) - spot_z)
        xy, heights = unique_samples(
            np.concatenate([xy, spot_xy]),
            np.concatenate([heights, spot_z]),
            **conflict_options,
        )
    minx, miny, maxx, maxy = projected_sheet.bounds
    west, south = (
        math.floor(minx / resolution) * resolution,
        math.floor(miny / resolution) * resolution,
    )
    east, north = (
        math.ceil(maxx / resolution) * resolution,
        math.ceil(maxy / resolution) * resolution,
    )
    width, height = (
        round((east - west) / resolution),
        round((north - south) / resolution),
    )
    if width * height > 25_000_000:
        raise ValueError("Pilot limited to 25 million cells; use a coarser resolution")
    grid_transform = from_origin(west, north, resolution, resolution)
    x = west + (np.arange(width) + 0.5) * resolution
    y = north - (np.arange(height) + 0.5) * resolution
    xx, yy = np.meshgrid(x, y)
    grid_xy = np.column_stack([xx.ravel(), yy.ravel()])
    inside = geometry_mask(
        [mapping(projected_sheet)], (height, width), grid_transform, invert=True
    )
    engine_report = None
    if interpolation == "imagico":
        if spots and validate_points:
            validation, _, _ = build_imagico_grid(
                interpolation_contours,
                [],
                unknown_contours,
                (height, width),
                grid_transform,
                progress=lambda stage: print(
                    f"Imagico validation (points held out): {stage}", flush=True
                ),
            )
            held_out = error_stats(
                sample_grid(validation.elevation, grid_transform, spot_xy) - spot_z
            )
            held_out.update(method="imagico_python", final_run_includes_all_points=True)
            validation.elevation[~inside] = np.nan
            write_raster(
                output_dir / "dem_contours_only.tif",
                validation.elevation,
                grid_transform,
                epsg,
                "Validation only: point heights held out",
                method="imagico_python_holdout",
            )
            (output_dir / "imagico_validation_report.json").write_text(
                json.dumps(
                    dict(
                        validation.report,
                        unknown_height_estimates=validation.unknown_heights,
                        point_errors=held_out,
                    ),
                    indent=2,
                )
                + "\n"
            )
            del validation
        result, known_grid, unknown_grid = build_imagico_grid(
            interpolation_contours,
            spots,
            unknown_contours,
            (height, width),
            grid_transform,
            progress=lambda stage: print(
                f"Imagico final (all points included): {stage}", flush=True
            ),
        )
        dem = result.elevation
        engine_report = dict(
            result.report, unknown_height_estimates=result.unknown_heights
        )
        np.savez_compressed(
            output_dir / "interpolation_inputs.npz",
            known_heights=known_grid,
            unknown_ids=unknown_grid,
        )
        (output_dir / "imagico_report.json").write_text(
            json.dumps(engine_report, indent=2) + "\n"
        )
        estimates = {r["unknown_id"]: r for r in result.unknown_heights}
        final_contours = []
        for line, props in recovered_contours:
            if "unknown_height_id" in props:
                estimate = estimates[props["unknown_height_id"]]
                props = dict(
                    props,
                    VALUE=estimate["estimated_elevation_m"],
                    height_status="estimated",
                    height_estimation=dict(
                        estimate, method="imagico_initial_surface_mean"
                    ),
                )
                props["height_recovery"] = dict(
                    props["height_recovery"],
                    recovered_elevation_m=estimate["estimated_elevation_m"],
                    method="imagico_unknown_height",
                )
            final_contours.append((line, props))
        recovered_contours = final_contours
        write_vectors(
            output_dir / "recovered_contours.geojsonl", recovered_contours, unproject
        )
        height_handling["estimated_contours"] = len(estimates)
        height_handling["estimates"] = result.unknown_heights
        classification_path = output_dir / "classification.json"
        classification = json.loads(classification_path.read_text())
        classification["height_recovery"] = height_handling
        classification_path.write_text(json.dumps(classification, indent=2) + "\n")

        def model(coordinates):
            return sample_grid(dem, grid_transform, coordinates)
    else:
        model = LinearNDInterpolator(xy, heights) if spots else contour_model
        dem = model(grid_xy).reshape(height, width)
    dem[~inside] = np.nan
    if not np.any(np.isfinite(dem)):
        raise ValueError("Interpolation produced no supported cells")
    distances = cKDTree(xy).query(grid_xy)[0].reshape(height, width)
    distances[~inside] = np.nan
    write_raster(
        output_dir / "dem.tif",
        dem,
        grid_transform,
        epsg,
        "Elevation above source vertical datum",
        method="imagico_python" if interpolation == "imagico" else "linear_TIN_pilot",
    )
    write_raster(
        output_dir / "constraint_distance.tif",
        distances,
        grid_transform,
        epsg,
        "Distance to nearest elevation sample",
    )
    water_decisions = []
    if water:
        polygons = [
            transform(project, shape(f["geometry"])) for f in read_features(water)
        ]
        flattened, water_decisions = flatten_water(
            dem, polygons, model, grid_transform, projected_sheet, resolution
        )
        write_raster(
            output_dir / "dem_water_flattened.tif",
            flattened,
            grid_transform,
            epsg,
            "Elevation; eligible water polygons flattened",
        )
    # This tests reproduction between sampled vertices, not independent terrain accuracy.
    check_xy, check_z = [], []
    for line, properties in interpolation_contours:
        dense = segmentize(line, resolution / 2)
        coords = np.asarray(dense.coords)[:, :2]
        midpoints = (coords[1:] + coords[:-1]) / 2
        check_xy.extend(midpoints)
        check_z.extend([properties["VALUE"]] * len(midpoints))
    residuals = error_stats(model(np.asarray(check_xy)) - np.asarray(check_z))
    report = {
        "sheet": sheet_id,
        "bounds_wgs84": list(sheet.bounds),
        "crs": f"EPSG:{epsg}",
        "resolution_m": resolution,
        "shape": [height, width],
        "method": "imagico_python"
        if interpolation == "imagico"
        else "unconstrained_linear_TIN",
        "interpolation_engine": engine_report,
        "source_files": [
            {"path": str(p), "sha256": hashlib.sha256(Path(p).read_bytes()).hexdigest()}
            for p in inputs
        ],
        "auxiliary_files": [
            {"path": str(p), "sha256": hashlib.sha256(Path(p).read_bytes()).hexdigest()}
            for p in (index, water, spot_review)
            if p is not None
        ],
        "counts": counts,
        "symbol_type_counts": dict(Counter(p["symbol_type"] for _, p in spots)),
        "point_layers": {
            "all": "spots.geojsonl",
            "triangles": "triangle_points.geojsonl",
        },
        "contour_elevation_counts": levels,
        "spot_radius_threshold_m": spot_radius,
        "spot_roundness_threshold": 0.98,
        "triangle_symbol_criteria": {
            "nominal_side_m": triangle_side,
            "side_range_m": [0.8 * triangle_side, 1.2 * triangle_side],
            "minimum_side_ratio": 0.8,
            "maximum_edge_error_m": triangle_side / 70.0,
            "maximum_area_error_fraction": 0.03,
        },
        "edge_sliver_tolerance_m": edge_sliver_m,
        "form_line_interval_m": form_line_interval_audit.get("interval_m")
        if form_line_interval == "auto"
        else form_line_interval,
        "form_line_detection": form_line_detection,
        "contour_contact_handling": contact_handling,
        "height_recovery": height_handling,
        "excluded_contours": "excluded_contours.geojsonl",
        "seam_reconciled": False,
        "constraint_samples": len(xy),
        "sheet_cells": int(inside.sum()),
        "valid_cells": int(np.isfinite(dem).sum()),
        "coverage_fraction": float(np.isfinite(dem).sum() / inside.sum()),
        "elevation_min_m": float(np.nanmin(dem)),
        "elevation_max_m": float(np.nanmax(dem)),
        "nearest_sample_distance_p95_m": float(np.nanpercentile(distances[inside], 95)),
        "contour_midpoint_residuals": residuals,
        "spot_holdout_contours_only": held_out,
        "point_holdout_requested": validate_points,
        "water": water_decisions,
        "limitations": [
            "30 m is grid spacing, not a claim of terrain accuracy.",
            "Imagico Python uses raster constraints and Gaussian filters; it is an adaptation, not a bit-for-bit CImg reproduction."
            if interpolation == "imagico"
            else "Unconstrained TIN may terrace hilltops, cross contour segments, and bridge data gaps.",
            "Circular and triangular spot symbols are classified heuristically; inspect spots.geojsonl.",
            "Source vertical datum and accuracy have not been independently verified.",
            "Directional interpolation can extend beyond the vector sample hull; no-hit cells remain masked. Sheet boundaries need later seam handling."
            if interpolation == "imagico"
            else "No extrapolation beyond sample convex hull; sheet boundaries need later seam handling.",
            "Outside edge sliver removal does not reconcile source height jumps across sheet seams.",
            "Optional form-line detection uses a calibrated junction heuristic; it is not an exhaustive subtype classifier.",
            "Small cartographic contacts relax local contour geometry; their true ground separation is not recovered.",
            "Height recovery uses vector heuristics; original heights and each correction are retained in the exported audit.",
            "Water levels are estimated from the interpolated surface, not measured; variable or incomplete shores are skipped.",
        ],
    }
    (output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: report[k]
                for k in [
                    "crs",
                    "shape",
                    "coverage_fraction",
                    "elevation_min_m",
                    "elevation_max_m",
                    "spot_holdout_contours_only",
                ]
            },
            indent=2,
        )
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "inputs", nargs="+", type=Path, help="iomaps GeoJSONL contour extracts (WGS84)"
    )
    parser.add_argument("--sheet", required=True)
    parser.add_argument(
        "--index", type=Path, default=Path(__file__).with_name("index_50k.geojson")
    )
    parser.add_argument(
        "--output-dir", required=True, type=Path, help="New directory for this run"
    )
    parser.add_argument("--resolution", type=float, default=30)
    parser.add_argument("--skip-holdout", action="store_true", help="Production: interpolate once with all points; omit the calibration-only holdout pass")
    parser.add_argument(
        "--interpolation",
        choices=["imagico", "linear"],
        default="imagico",
        help="Imagico Python stages (default), or the historical linear TIN baseline",
    )
    parser.add_argument(
        "--spot-radius",
        type=float,
        default=13,
        help="Maximum circular symbol radius in metres; 0 disables",
    )
    parser.add_argument(
        "--water", type=Path, help="Optional iomaps CartoDEM waterbody GeoJSONL"
    )
    parser.add_argument(
        "--spot-review",
        type=Path,
        help="Optional geometry-verified spot corrections JSON",
    )
    parser.add_argument(
        "--triangle-side",
        type=float,
        default=70,
        help="Nominal triangular symbol side in metres (20%% tolerance); 0 disables",
    )
    parser.add_argument(
        "--edge-sliver-m",
        type=float,
        default=2.0,
        help="Exclude mostly-outside lines whose retained portion lies within this edge band; 0 disables",
    )
    parser.add_argument(
        "--contour-contact-radius-m",
        type=float,
        help="Relax local continuing-contour contacts within this radius; default: half grid spacing; 0 disables",
    )
    parser.add_argument(
        "--form-line-interval",
        type=parse_form_line_interval,
        default="auto",
        help="auto (default): infer one interval per sheet; none: disable; positive metres: supply a known interval",
    )
    parser.add_argument(
        "--height-recovery",
        choices=["auto", "none"],
        default="auto",
        help="Recover supported closed-loop attribute errors; none disables",
    )
    parser.add_argument(
        "--height-repair-max-length-m",
        type=float,
        default=10.0,
        help="Maximum total fragment length for the small-piece fallback; 0 disables that fallback",
    )
    parser.add_argument(
        "--height-repair-max-fraction",
        type=float,
        default=0.10,
        help="Maximum fraction of a loop changed by the small-piece fallback (default 0.10)",
    )
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("Output directory is not empty; choose a new run directory")
    try:
        build_dem(
            args.inputs,
            args.index,
            args.sheet,
            args.output_dir,
            args.resolution,
            args.water,
            args.spot_radius,
            args.spot_review,
            args.triangle_side,
            edge_sliver_m=args.edge_sliver_m,
            form_line_interval=args.form_line_interval,
            contour_contact_radius_m=args.contour_contact_radius_m,
            height_recovery=args.height_recovery,
            height_repair_max_length_m=args.height_repair_max_length_m,
            height_repair_max_fraction=args.height_repair_max_fraction,
            interpolation=args.interpolation,
            validate_points=not args.skip_holdout,
        )
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
