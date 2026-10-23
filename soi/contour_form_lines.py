"""Infer one regular contour interval per sheet for the calibrated SOI form-line junction pattern.

This is a geometry heuristic, not a general form-line subtype classifier.
Uncertain branches remain DEM constraints, including any height conflict.
"""

import math
from collections import Counter, defaultdict
from itertools import pairwise

from shapely import STRtree, line_merge
from shapely.geometry import LineString, MultiLineString, Point


def parse_form_line_interval(value):
    if value is None or value == "none":
        return None
    if value == "auto":
        return "auto"
    try:
        value = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(
            "Form-line interval must be auto, none, or positive metres"
        ) from exc
    if not math.isfinite(value) or value <= 0:
        raise ValueError("Form-line contour interval must be positive and finite")
    return value


def contour_chains(records):
    """Join only exact same-height endpoints; preserve every source member."""
    groups = defaultdict(list)
    for number, (_, _, metadata) in enumerate(records):
        groups[metadata["VALUE"]].append(number)
    chains = []
    for elevation, indices in groups.items():
        lines = [records[i][0] for i in indices]
        tree = STRtree(lines)
        merged = line_merge(MultiLineString(lines))
        for line in [merged] if merged.geom_type == "LineString" else merged.geoms:
            members = sorted(
                indices[int(i)] for i in tree.query(line, predicate="covers")
            )
            metric_coordinates = {
                xy: projected_xy
                for i in members
                for xy, projected_xy in zip(records[i][0].coords, records[i][1].coords)
            }
            projected = LineString([metric_coordinates[xy] for xy in line.coords])
            chains.append((line, projected, elevation, members))
    return chains


def infer_sheet_interval(chains, branches, sheet=None):
    """Infer the dominant regular series using contour length inside the sheet.

    Point symbols and suspect branches cannot vote. Small off-series loops and
    attribute fragments remain constraints but need not define the interval.
    """
    lengths = defaultdict(float)
    for number, (_, line, height, _) in enumerate(chains):
        if number not in branches:
            length = line.length if sheet is None else line.intersection(sheet).length
            if length > 0:
                lengths[height] += length
    minimum = 0.1 * max(lengths.values(), default=0)
    strong = {z: length for z, length in lengths.items() if length >= minimum}
    evidence = {
        "scope": "selected_sheet",
        "assumption": "one_regular_interval_per_sheet",
        "minimum_reference_length_m": minimum,
        "reference_lengths_m": [
            {"elevation_m": z, "length_m": length, "used_for_interval": z in strong}
            for z, length in sorted(lengths.items())
        ],
        "status": "insufficient_levels",
        "interval_m": None,
        "reference_elevation_m": None,
    }
    if len(strong) < 3:
        return evidence
    votes = Counter()
    for lower, upper in pairwise(sorted(strong)):
        gap = round(upper - lower, 8)
        if gap > 0:
            votes[gap] += min(strong[lower], strong[upper])
    if not votes:
        return evidence
    interval, weight = votes.most_common(1)[0]
    anchor = max(strong, key=strong.get)
    aligned = {
        z: length
        for z, length in strong.items()
        if math.isclose(
            (z - anchor) / interval,
            round((z - anchor) / interval),
            abs_tol=1e-8,
            rel_tol=0,
        )
    }
    support = weight / sum(votes.values())
    alignment = sum(aligned.values()) / sum(strong.values())
    supported = support > 0.5 and alignment >= 0.9 and len(aligned) >= 3
    evidence.update(
        gap_support_fraction=support,
        reference_lattice_support_fraction=alignment,
        reference_elevation_m=anchor,
        status="supported" if supported else "mixed_spacing",
        interval_m=interval if supported else None,
    )
    return evidence


def identify_form_line_candidates(
    records, interval_m, *, audit=None, sheet=None, interval_audit=None
):
    """Classify open branches ending at exact vertices of continuing contours.

    In auto mode, infer one interval from the selected sheet's linework. A
    branch must be half that interval below its parent at every junction. Known
    numeric intervals remain available; neither mode assumes elevations start at 0.
    """
    mode = parse_form_line_interval(interval_m)
    if mode is None or not records:
        if interval_audit is not None:
            interval_audit.update(
                status="disabled" if mode is None else "no_contours", interval_m=None
            )
        return {}
    chains = contour_chains(records)
    tree = STRtree([c[0] for c in chains])
    contacts = defaultdict(list)
    for number, (line, projected, height, _) in enumerate(chains):
        if line.is_closed:
            continue
        for endpoint, metric_xy in zip(
            (line.coords[0], line.coords[-1]),
            (projected.coords[0], projected.coords[-1]),
        ):
            for other in tree.query(Point(endpoint), predicate="intersects"):
                other = int(other)
                target, _, target_height, _ = chains[other]
                if other == number or target_height == height:
                    continue
                if not target.is_closed and endpoint in (
                    target.coords[0],
                    target.coords[-1],
                ):
                    continue
                if endpoint not in target.coords:
                    continue  # Do not snap or invent a vertex at an unsampled crossing.
                contacts[number].append((other, endpoint, metric_xy))
    inference = (
        infer_sheet_interval(chains, contacts, sheet)
        if mode == "auto"
        else {"scope": "selected_sheet", "status": "supplied", "interval_m": mode}
    )
    if interval_audit is not None:
        interval_audit.update(inference)
    interval = inference["interval_m"]
    candidates = {}
    for number, junctions in contacts.items():
        _, _, height, members = chains[number]
        evidence = {
            "classification_method": "half_interval_branch_to_continuing_chain",
            "interval_source": "sheet_contour_spacing"
            if mode == "auto"
            else "supplied",
            "chain_source_geometry_sha256": [
                records[i][2]["source_geometry_sha256"] for i in members
            ],
            "junctions": [],
        }
        for other, xy, metric_xy in junctions:
            target_height = chains[other][2]
            evidence["junctions"].append(
                {
                    "coordinates": list(xy),
                    "elevation_m": target_height,
                    "source_ids": [
                        records[i][2].get("source_id") for i in chains[other][3]
                    ],
                    "source_geometry_sha256": [
                        records[i][2]["source_geometry_sha256"]
                        for i in chains[other][3]
                    ],
                    "regular_contour_interval_m": interval,
                }
            )
        if interval is None:
            decision = "retain_uncertain_interval"
        elif mode == "auto" and any(
            not math.isclose(
                (j["elevation_m"] - inference["reference_elevation_m"]) / interval,
                round(
                    (j["elevation_m"] - inference["reference_elevation_m"]) / interval
                ),
                rel_tol=0,
                abs_tol=1e-8,
            )
            for j in evidence["junctions"]
        ):
            decision = "retain_off_series_parent"
        elif all(
            math.isclose(
                j["elevation_m"] - height, interval / 2, abs_tol=1e-8, rel_tol=0
            )
            for j in evidence["junctions"]
        ):
            decision = "exclude_form_line_candidate"
        else:
            decision = "retain_non_half_interval_join"
        evidence["decision"] = decision
        evidence["regular_contour_interval_m"] = interval
        evidence["sheet_interval_status"] = inference["status"]
        for member in members:
            if decision == "exclude_form_line_candidate":
                candidates[member] = evidence
            if audit is not None:
                _, projected, metadata = records[member]
                audit.append((projected, dict(metadata, form_line_evidence=evidence)))
    return candidates
