# /// script
# requires-python = ">=3.12"
# dependencies = ["shapely>=2.1", "pytest>=8"]
# ///
import hashlib

import pytest
from shapely.geometry import LineString, box

from contour_form_lines import identify_form_line_candidates, parse_form_line_interval


def records(lines):
    return [
        (
            line,
            line,
            {
                "VALUE": height,
                "source_id": str(i),
                "source_geometry_sha256": hashlib.sha256(line.wkb).hexdigest(),
            },
        )
        for i, (line, height) in enumerate(lines)
    ]


def landscape(interval, origin=100, x=0, branch_step=0.5):
    lines = [
        (
            LineString([(x + i * 180, -4000), (x + i * 180, 0), (x + i * 180, 4000)]),
            origin + i * interval,
        )
        for i in range(-2, 3)
    ]
    lines.append((LineString([(x, 0), (x - 90, 100)]), origin - branch_step * interval))
    return lines


@pytest.mark.parametrize("interval", [5, 10, 20, 40])
@pytest.mark.parametrize("offset", [0, 3.7, -143.25])
def test_sheet_spacing_adapts_and_does_not_depend_on_height_origin(interval, offset):
    data = records(landscape(interval, 100 + offset))
    audit, inference = [], {}
    found = identify_form_line_candidates(
        data, "auto", audit=audit, interval_audit=inference
    )
    assert set(found) == {5}
    assert found[5]["regular_contour_interval_m"] == pytest.approx(interval)
    assert inference["status"] == "supported"
    assert inference["scope"] == "selected_sheet"
    assert audit[0][1]["VALUE"] == pytest.approx(100 + offset - interval / 2)


def test_finer_sheet_spacing_keeps_full_interval_branches():
    # A solid 90 m contour touching 100 m must not be mistaken for a form line
    # because most of the sheet uses a coarser interval elsewhere.
    data = records(landscape(10, branch_step=1) + landscape(20, x=10000))
    found = identify_form_line_candidates(data, "auto")
    assert found == {}


def test_mixed_spacing_and_sparse_evidence_keep_conflicting_lines():
    sparse = [landscape(20)[2], landscape(20)[-1]]
    audit = []
    assert identify_form_line_candidates(records(sparse), "auto", audit=audit) == {}
    assert audit[0][1]["form_line_evidence"]["decision"] == "retain_uncertain_interval"
    mixed = landscape(20)
    mixed.insert(0, (LineString([(90, -4000), (90, 4000)]), 110))
    assert identify_form_line_candidates(records(mixed), "auto") == {}


def test_minor_off_series_lines_do_not_override_the_sheet_interval():
    lines = landscape(20)
    lines.extend(
        [
            (LineString([(60, -700), (60, 700)]), 110),
            (LineString([(-60, -700), (-60, 700)]), 90),
        ]
    )
    assert set(identify_form_line_candidates(records(lines), "auto")) == {5}


def test_all_junctions_use_the_same_sheet_interval():
    lines = landscape(20)[:-1] + landscape(40, origin=110, x=10000)[:-1]
    lines.append((LineString([(0, 0), (5000, 100), (10000, 0)]), 90))
    audit = []
    assert identify_form_line_candidates(records(lines), "auto", audit=audit) == {}
    evidence = audit[0][1]["form_line_evidence"]
    assert evidence["decision"] == "retain_non_half_interval_join"
    assert {j["regular_contour_interval_m"] for j in evidence["junctions"]} == {10}


def test_distant_reference_contours_support_sparse_junction_but_other_sheets_do_not():
    lines = landscape(20)[:-1] + [
        (LineString([(10000, -4000), (10000, 0), (10000, 4000)]), 100),
        (LineString([(10000, 0), (9900, 100)]), 90),
    ]
    assert set(identify_form_line_candidates(records(lines), "auto")) == {6}
    # An extract can contain uncut lines from outside the requested sheet.
    assert (
        identify_form_line_candidates(
            records(lines), "auto", sheet=box(9000, -5000, 11000, 5000)
        )
        == {}
    )


def test_exact_fragment_joins_and_direction_do_not_change_classification():
    lines = landscape(20)
    fragmented = []
    for line, height in lines:
        points = list(line.coords)
        if len(points) == 3:
            fragmented.extend(
                [
                    (LineString(points[:2]), height),
                    (LineString(points[1:][::-1]), height),
                ]
            )
        else:
            mid = (-45, 50)
            fragmented.extend(
                [
                    (LineString([points[0], mid]), height),
                    (LineString([points[-1], mid]), height),
                ]
            )
    found = identify_form_line_candidates(records(fragmented), "auto")
    assert set(found) == {10, 11}
    assert all(e["regular_contour_interval_m"] == 20 for e in found.values())
    assert len(found[10]["chain_source_geometry_sha256"]) == 2


def test_closed_parent_start_vertex_and_known_offset_interval():
    ring = LineString([(0, 0), (200, 0), (200, 200), (0, 200), (0, 0)])
    branch = LineString([(0, 0), (-100, -100)])
    data = records([(ring, 103.7), (branch, 93.7)])
    assert set(identify_form_line_candidates(data, 20)) == {1}
    assert identify_form_line_candidates(data, 10) == {}
    assert identify_form_line_candidates(data, None) == {}


def test_endpoint_joins_and_near_misses_are_not_form_line_evidence():
    lines = landscape(20)[:-1]
    lines.extend(
        [
            (LineString([(0, 4000), (100, 4100)]), 90),
            (LineString([(0.001, 0), (-90, 100)]), 90),
        ]
    )
    assert identify_form_line_candidates(records(lines), "auto") == {}


@pytest.mark.parametrize("value", [0, -10, "nan", "inf", "bogus"])
def test_invalid_intervals_are_rejected(value):
    with pytest.raises(ValueError, match="interval"):
        parse_form_line_interval(value)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
