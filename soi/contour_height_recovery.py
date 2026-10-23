"""Propose contour-height repairs from closed-loop topology and sheet spacing.

This pass consumes classified contours, never raster pixels or
spot-height symbols. It returns a derived copy and an audit; inputs stay intact.
"""

import math

from shapely import STRtree, line_merge
from shapely.geometry import MultiLineString, Point, Polygon

from contour_form_lines import contour_chains


def validate_recovery_limits(max_fragment_length_m, max_fragment_fraction):
    if not math.isfinite(max_fragment_length_m) or max_fragment_length_m < 0:
        raise ValueError("Height repair length must be finite and nonnegative")
    if not math.isfinite(max_fragment_fraction) or not 0 <= max_fragment_fraction < 0.5:
        raise ValueError("Height repair fraction must be finite and in [0, 0.5)")


def recover_contour_heights(
    contours,
    sheet,
    sheet_interval,
    *,
    boundary_margin_m=30,
    max_fragment_length_m=10.0,
    max_fragment_fraction=0.10,
):
    """Recover two-height loops from nesting evidence or a size-limited fragment.

    Two open, same-height chains must join each other at both ends and form a
    simple closed contour. Its nearest enclosing regular contour and its
    immediate regular inner contours must bracket one existing height by
    exactly one inferred sheet interval. Crossings, contacts, conflicting inner
    loops and proximity to the sheet edge prevent recovery.
    When references are missing or incomplete, a small-fragment fallback accepts
    at most 10 m and 10% of the whole loop by default. The replacement is the
    existing regular height joined at both ends. Length is measured across all
    connected same-height fragments, never separately per source feature.

    Returns (derived_contours, cycle_audit). This is a single pass: a proposed
    correction never becomes evidence for another correction in the same run.
    """
    if not math.isfinite(boundary_margin_m) or boundary_margin_m < 0:
        raise ValueError("Boundary margin must be finite and nonnegative")
    validate_recovery_limits(max_fragment_length_m, max_fragment_fraction)
    if not contours:
        return contours, []
    chains = contour_chains([(g, g, p) for g, p in contours])
    tree = STRtree([row[0] for row in chains])
    cycles, seen = [], set()
    for number, (line, _, height, members) in enumerate(chains):
        if line.is_closed or not line.is_simple:
            continue
        others = []
        for xy in (line.coords[0], line.coords[-1]):
            neighbours = [
                int(i)
                for i in tree.query(Point(xy), predicate="intersects")
                if int(i) != number
            ]
            if len(neighbours) != 1:
                break
            other = neighbours[0]
            other_line = chains[other][0]
            if xy not in (other_line.coords[0], other_line.coords[-1]):
                break
            others.append(other)
        if len(others) != 2 or others[0] != others[1]:
            continue
        pair = tuple(sorted((number, others[0])))
        if pair in seen or height == chains[others[0]][2]:
            continue
        seen.add(pair)
        whole = line_merge(MultiLineString([chains[i][0] for i in pair]))
        if whole.geom_type != "LineString" or not whole.is_ring:
            continue
        polygon = Polygon(whole)
        if not polygon.is_valid or polygon.area == 0:
            continue
        cycles.append((whole, polygon, pair))

    interval = sheet_interval.get("interval_m")
    anchor = sheet_interval.get("reference_elevation_m")
    supported = (
        sheet_interval.get("status") == "supported"
        and interval is not None
        and math.isfinite(interval)
        and interval > 0
        and anchor is not None
        and math.isfinite(anchor)
    )

    def aligned(height):
        return supported and math.isclose(
            (height - anchor) / interval,
            round((height - anchor) / interval),
            abs_tol=1e-8,
            rel_tol=0,
        )

    rings = [
        (Polygon(line), number, height)
        for number, (line, _, height, _) in enumerate(chains)
        if line.is_ring and aligned(height) and Polygon(line).is_valid
    ]
    ring_tree = STRtree([row[0] for row in rings])

    def source_ids(numbers):
        return sorted(
            {
                contours[m][1].get("source_id", contours[m][1].get("id", str(m)))
                for n in numbers
                for m in chains[n][3]
            }
        )

    def describe_reference(row):
        return {"elevation_m": row[2], "source_ids": source_ids([row[1]])}

    audit, changes = [], {}
    for whole, polygon, pair in cycles:
        evidence = {
            "recovery_id": f"H{len(audit) + 1:04d}",
            "method": "closed_loop_between_regular_contours",
            "original_elevations_m": sorted(chains[n][2] for n in pair),
            "proposed_elevation_m": None,
            "regular_contour_interval_m": interval,
            "source_ids": source_ids(pair),
            "boundary_distance_m": whole.distance(sheet.boundary),
            "cycle_length_m": whole.length,
            "decision": "retain_uncertain_interval",
        }
        audit.append((whole, evidence))
        if not supported:
            continue
        if (
            not sheet.contains(polygon)
            or evidence["boundary_distance_m"] <= boundary_margin_m
        ):
            evidence["decision"] = "retain_near_sheet_boundary"
            continue
        contacts = [
            int(i)
            for i in tree.query(whole, predicate="intersects")
            if int(i) not in pair
        ]
        if contacts:
            evidence["decision"] = "retain_contact_or_crossing"
            continue
        parents = [rings[int(i)] for i in ring_tree.query(polygon, predicate="within")]
        parent = min(parents, key=lambda row: row[0].area) if parents else None
        enclosed = [
            rings[int(i)] for i in ring_tree.query(polygon, predicate="contains")
        ]
        children = []
        for child in sorted(enclosed, key=lambda row: -row[0].area):
            if not any(row[0].contains(child[0]) for row in children):
                children.append(child)
        evidence["outer_reference"] = describe_reference(parent) if parent else None
        evidence["inner_references"] = [describe_reference(row) for row in children]
        if parent is None or not children:
            evidence["decision"] = "retain_missing_reference"
            continue
        if any(
            int(i) != row[1]
            for row in [parent, *children]
            for i in tree.query(chains[row[1]][0], predicate="intersects")
        ):
            evidence["decision"] = "retain_reference_contact_or_crossing"
            continue
        # Another unresolved cycle between the candidate and its references is
        # a missing immediate reference, not permission to jump over a contour.
        if any(
            other_pair != pair
            and (
                (other.contains(polygon) and parent[0].contains(other))
                or (
                    polygon.contains(other)
                    and not any(c[0].contains(other) for c in children)
                )
            )
            for _, other, other_pair in cycles
        ):
            evidence["decision"] = "retain_unresolved_neighbour"
            continue
        # Open linework between this loop and its nearest inner closed contours
        # means that nesting alone does not describe all the available evidence.
        if any(
            int(i) not in pair
            and not chains[int(i)][0].is_ring
            and not any(c[0].covers(chains[int(i)][0]) for c in children)
            for i in tree.query(polygon, predicate="contains")
        ):
            evidence["decision"] = "retain_open_inner_reference"
            continue
        # The same completeness check applies outside the candidate. Other
        # closed regular loops can protect unrelated linework in sibling hills.
        siblings = [
            rings[int(i)][0]
            for i in ring_tree.query(parent[0], predicate="contains")
            if rings[int(i)][0].disjoint(polygon)
        ]
        if any(
            int(i) not in pair
            and not chains[int(i)][0].is_ring
            and not polygon.covers(chains[int(i)][0])
            and not any(sibling.covers(chains[int(i)][0]) for sibling in siblings)
            for i in tree.query(parent[0], predicate="contains")
        ):
            evidence["decision"] = "retain_open_outer_reference"
            continue
        inner_heights = {row[2] for row in children}
        if len(inner_heights) != 1:
            evidence["decision"] = "retain_ambiguous_nesting"
            continue
        outer, inner = parent[2], next(iter(inner_heights))
        if not math.isclose(abs(inner - outer), 2 * interval, abs_tol=1e-8, rel_tol=0):
            evidence["decision"] = "retain_ambiguous_nesting"
            continue
        proposed = (inner + outer) / 2
        evidence["proposed_elevation_m"] = proposed
        matched = [
            chains[n][2]
            for n in pair
            if math.isclose(chains[n][2], proposed, abs_tol=1e-8, rel_tol=0)
        ]
        if not matched:
            evidence["decision"] = "retain_neither_height_supported"
            continue
        # Keep the exact existing numeric value so harmless floating-point
        # arithmetic cannot introduce another tiny mismatch at joined endpoints.
        proposed = matched[0]
        evidence["proposed_elevation_m"] = proposed
        # Off-series closed loops do not set the regular interval. Direct inner
        # off-series contours must still lie in the same direction of relief.
        off_series = [
            chains[int(i)][2]
            for i in tree.query(polygon, predicate="contains")
            if int(i) not in pair
            and not aligned(chains[int(i)][2])
            and not any(c[0].covers(chains[int(i)][0]) for c in children)
        ]
        evidence["direct_inner_off_series_elevations_m"] = sorted(set(off_series))
        direction = 1 if inner > outer else -1
        outer_off_series = [
            height
            for line, _, height, _ in chains
            if line.is_ring
            and not aligned(height)
            and parent[0].contains(line)
            and Polygon(line).contains(polygon)
        ]
        evidence["direct_outer_off_series_elevations_m"] = sorted(set(outer_off_series))
        if any(
            not 0 < direction * (z - proposed) < interval for z in off_series
        ) or any(
            not 0 < direction * (proposed - z) < interval for z in outer_off_series
        ):
            evidence["decision"] = "retain_contradictory_off_series_reference"
            continue
        evidence["decision"] = "supported_recovery"
        changed = [m for n in pair if chains[n][2] != proposed for m in chains[n][3]]
        evidence["changed_source_ids"] = sorted(
            {contours[m][1].get("source_id", str(m)) for m in changed}
        )
        evidence["changed_length_m"] = sum(contours[m][0].length for m in changed)
        for member in changed:
            changes[member] = evidence
    # Apply this deliberately limited fallback only after the stronger nesting
    # checks. Known contradictory references, boundary cases and contacts never
    # become eligible merely because a fragment is short.
    for (_, evidence), (whole, polygon, pair) in zip(audit, cycles):
        if evidence["decision"] not in {
            "retain_missing_reference",
            "retain_open_inner_reference",
            "retain_open_outer_reference",
        }:
            continue
        small, rest = sorted(pair, key=lambda n: chains[n][0].length)
        length = chains[small][0].length
        fraction = length / whole.length
        fragment = {
            "length_m": length,
            "fraction_of_loop": fraction,
            "maximum_length_m": max_fragment_length_m,
            "maximum_fraction_of_loop": max_fragment_fraction,
            "replacement_elevation_m": chains[rest][2],
            "decision": "retain_size_limit",
        }
        evidence["fragment_recovery"] = fragment
        if max_fragment_length_m == 0 or max_fragment_fraction == 0:
            fragment["decision"] = "disabled"
            continue
        if (
            length > max_fragment_length_m + 1e-8
            or fraction > max_fragment_fraction + 1e-10
        ):
            continue
        if not aligned(chains[rest][2]):
            fragment["decision"] = "retain_off_series_replacement"
            continue
        replacement = chains[rest][2]
        outer_reference = evidence.get("outer_reference")
        inner_references = evidence.get("inner_references", [])
        known = ([outer_reference] if outer_reference else []) + inner_references
        inner_values = {r["elevation_m"] for r in inner_references}
        if (
            any(
                not math.isclose(
                    abs(r["elevation_m"] - replacement),
                    interval,
                    abs_tol=1e-8,
                    rel_tol=0,
                )
                for r in known
            )
            or len(inner_values) > 1
            or (
                outer_reference
                and inner_values
                and not math.isclose(
                    replacement - outer_reference["elevation_m"],
                    next(iter(inner_values)) - replacement,
                    abs_tol=1e-8,
                    rel_tol=0,
                )
            )
        ):
            fragment["decision"] = "retain_contradictory_reference"
            continue
        evidence["nesting_decision"] = evidence["decision"]
        evidence.update(
            decision="supported_recovery",
            method="small_fragment_between_same_height_endpoints",
            proposed_elevation_m=chains[rest][2],
            changed_source_ids=source_ids([small]),
            changed_length_m=length,
        )
        fragment["decision"] = "recover_within_size_limits"
        for member in chains[small][3]:
            changes[member] = evidence
    derived = []
    for number, (line, props) in enumerate(contours):
        if number not in changes:
            derived.append((line, props))
            continue
        evidence = changes[number]
        derived.append(
            (
                line,
                dict(
                    props,
                    VALUE=evidence["proposed_elevation_m"],
                    height_recovery={
                        "recovery_id": evidence["recovery_id"],
                        "original_elevation_m": props["VALUE"],
                        "recovered_elevation_m": evidence["proposed_elevation_m"],
                        "method": evidence["method"],
                    },
                ),
            )
        )
    return derived, audit
