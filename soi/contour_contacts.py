"""Relax small cartographic contacts while preserving both contour elevations."""

import math
from collections import defaultdict
from itertools import combinations, pairwise

from shapely.geometry import Point
from shapely.ops import unary_union

from contour_form_lines import contour_chains


def line_parts(geometry):
    if geometry.is_empty:
        return
    if geometry.geom_type == "LineString":
        yield geometry
    elif hasattr(geometry, "geoms"):
        for component in geometry.geoms:
            yield from line_parts(component)


def has_conflicting_vertices(contours):
    heights = {}
    for line, props in contours:
        elevation = props["VALUE"]
        for xy in line.coords:
            key = tuple(round(v, 3) for v in xy)
            previous = heights.setdefault(key, elevation)
            if abs(previous - elevation) > 0.01:
                return True
    return False


def relax_contour_contacts(contours, radius_m, sheet):
    """Return interpolation lines, contact audit, and relaxed source sections.

    Both heights must continue through a node. Accept isolated tangencies or
    paired crossings forming a lens no longer than twice the radius along either
    contour. Endpoint height changes, T-junctions, isolated crossings, overlaps,
    and contacts near the sheet boundary remain strict conflicts. No contour
    interval is assumed. Point heights are not inputs to this operation.
    """
    if not math.isfinite(radius_m) or radius_m < 0:
        raise ValueError("Contour contact radius must be nonnegative and finite")
    if radius_m == 0 or not contours:
        return contours, [], []
    if not has_conflicting_vertices(contours):
        return contours, [], []
    chains = contour_chains([(line, line, props) for line, props in contours])
    nodes = defaultdict(lambda: defaultdict(lambda: {"chains": set(), "rays": set()}))
    positions = {}
    for number, (line, _, height, _) in enumerate(chains):
        for a, b in pairwise(line.coords):
            ka, kb = tuple(round(v, 3) for v in a), tuple(round(v, 3) for v in b)
            if ka == kb:
                continue
            positions.setdefault(ka, a)
            positions.setdefault(kb, b)
            for here, there in ((ka, kb), (kb, ka)):
                nodes[here][height]["chains"].add(number)
                nodes[here][height]["rays"].add(there)
    candidates, by_pair = [], defaultdict(list)
    for xy, heights in nodes.items():
        if len(heights) != 2 or any(
            len(h["rays"]) != 2 or len(h["chains"]) != 1 for h in heights.values()
        ):
            continue
        point = Point(positions[xy])
        if point.distance(sheet.boundary) <= radius_m:
            continue
        rays = sorted(
            (math.atan2(p[1] - xy[1], p[0] - xy[0]) % (2 * math.pi), z)
            for z, data in heights.items()
            for p in data["rays"]
        )
        # Coincident rays imply a shared segment, not an isolated contact.
        angular_gaps = [
            (b[0] - a[0]) % (2 * math.pi) for a, b in zip(rays, rays[1:] + rays[:1])
        ]
        if min(angular_gaps) < 1e-8:
            continue
        order = [z for _, z in rays]
        changes = sum(a != b for a, b in zip(order, order[1:] + order[:1]))
        pair = tuple(sorted(next(iter(h["chains"])) for h in heights.values()))
        if any(not chains[n][0].is_simple for n in pair):
            continue
        candidate = {
            "point": point,
            "chains": pair,
            "kind": "tangent_contact" if changes == 2 else None,
        }
        by_pair[pair].append(len(candidates))
        candidates.append(candidate)
    for pair, indices in by_pair.items():
        for first, second in combinations(indices, 2):
            a, b = candidates[first], candidates[second]
            if a["point"].distance(b["point"]) > 2 * radius_m:
                continue
            distances = []
            for number in pair:
                line = chains[number][0]
                distance = abs(line.project(a["point"]) - line.project(b["point"]))
                if line.is_closed:
                    distance = min(distance, line.length - distance)
                distances.append(distance)
            if max(distances) <= 2 * radius_m:
                a["kind"] = b["kind"] = "short_contact_lens"
    active = {i for i, c in enumerate(candidates) if c["kind"] is not None}
    # A small closed contour must never disappear completely inside the relaxed
    # region. Preserve the entire connected contour if that would happen.
    while active:
        by_chain = defaultdict(list)
        for i in active:
            for number in candidates[i]["chains"]:
                by_chain[number].append(i)
        masks = {
            n: unary_union([candidates[i]["point"].buffer(radius_m) for i in ids])
            for n, ids in by_chain.items()
        }
        reject = {
            i
            for n, ids in by_chain.items()
            if chains[n][0].difference(masks[n]).length <= 0.001
            for i in ids
        }
        if not reject:
            break
        active -= reject
    if not active:
        return contours, [], []
    # Source exports remain intact. Only the interpolation copy loses the small
    # ambiguous sections; each surviving section keeps its original elevation.
    by_member = defaultdict(list)
    audit = []
    for i in sorted(active):
        candidate = candidates[i]
        members = sorted({m for n in candidate["chains"] for m in chains[n][3]})
        contact_id = f"K{len(audit) + 1:04d}"
        evidence = {
            "contact_id": contact_id,
            "classification": "cartographic_contact_candidate",
            "geometry_pattern": candidate["kind"],
            "handling": "relax_local_contour_geometry",
            "radius_m": radius_m,
            "elevations_m": sorted(chains[n][2] for n in candidate["chains"]),
            "source_geometry_sha256": sorted(
                {
                    contours[m][1].get("source_geometry_sha256")
                    for m in members
                    if contours[m][1].get("source_geometry_sha256")
                }
            ),
            "source_heights_changed": False,
        }
        audit.append((candidate["point"], evidence))
        region = candidate["point"].buffer(radius_m)
        for member in members:
            if contours[member][0].intersects(region):
                by_member[member].append((contact_id, region))
    retained, relaxed = [], []
    for number, (line, props) in enumerate(contours):
        if number not in by_member:
            retained.append((line, props))
            continue
        regions = by_member[number]
        mask = unary_union([region for _, region in regions])
        metadata = dict(props, contact_ids=[cid for cid, _ in regions])
        retained.extend((part, metadata) for part in line_parts(line.difference(mask)))
        relaxed.extend(
            (
                part,
                dict(
                    metadata,
                    constraint_status="relaxed_geometry",
                    reason="cartographic_contour_contact",
                ),
            )
            for part in line_parts(line.intersection(mask))
        )
    return retained, audit, relaxed
