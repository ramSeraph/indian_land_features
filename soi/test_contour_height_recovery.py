# /// script
# requires-python = ">=3.12"
# dependencies = ["shapely>=2.1", "pytest>=8"]
# ///
import copy

import pytest
from shapely.geometry import LineString, box

from contour_height_recovery import recover_contour_heights

SHEET = box(-1000, -1000, 1000, 1000)


def records(lines):
    return [
        (LineString(xy), {"VALUE": z, "source_id": f"source-{i}"})
        for i, (xy, z) in enumerate(lines)
    ]


def square(size):
    return [(-size, -size), (-size, size), (size, size), (size, -size), (-size, -size)]


def landscape(interval=20, origin=100, direction=1):
    target = origin + direction * interval
    inner = origin + direction * 2 * interval
    return records(
        [
            (square(100), origin),
            (square(60)[:4], target),
            ([(60, -60), (0, -60)], inner),
            ([(0, -60), (-60, -60)], inner),
            (square(20), inner),
        ]
    )


def inference(interval=20, origin=100):
    return {
        "status": "supported",
        "interval_m": interval,
        "reference_elevation_m": origin,
    }


@pytest.mark.parametrize("interval", [5, 10, 20, 40])
@pytest.mark.parametrize("origin", [100, -13.7])
@pytest.mark.parametrize("direction", [-1, 1])
def test_recovers_nested_hill_or_depression_at_inferred_interval(
    interval, origin, direction
):
    original = landscape(interval, origin, direction)
    snapshot = copy.deepcopy(original)
    derived, audit = recover_contour_heights(
        original, SHEET, inference(interval, origin)
    )
    assert len(audit) == 1
    assert audit[0][1]["decision"] == "supported_recovery"
    assert audit[0][1]["proposed_elevation_m"] == pytest.approx(
        origin + direction * interval
    )
    assert [i for i, (_, p) in enumerate(derived) if "height_recovery" in p] == [2, 3]
    assert original == snapshot
    assert all(a.wkb == b.wkb for (a, _), (b, _) in zip(original, derived))
    assert [p["VALUE"] for _, p in derived][1:4] == pytest.approx(
        [origin + direction * interval] * 3
    )
    for _, props in derived[2:4]:
        assert (
            props["height_recovery"]["original_elevation_m"] == original[2][1]["VALUE"]
        )
    assert recover_contour_heights(derived, SHEET, inference(interval, origin)) == (
        derived,
        [],
    )


def test_recovery_does_not_assume_the_shorter_fragment_is_wrong():
    original = landscape()
    original[1][1]["VALUE"] = 140
    original[2][1]["VALUE"] = original[3][1]["VALUE"] = 120
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert audit[0][1]["decision"] == "supported_recovery"
    assert [i for i, (_, p) in enumerate(derived) if "height_recovery" in p] == [1]


def test_fragment_order_and_direction_do_not_change_the_solution():
    original = landscape()
    reversed_parts = [(LineString(list(g.coords)[::-1]), p) for g, p in original[::-1]]
    a, _ = recover_contour_heights(original, SHEET, inference())
    b, _ = recover_contour_heights(reversed_parts, SHEET, inference())
    assert {p["source_id"]: p["VALUE"] for _, p in a} == {
        p["source_id"]: p["VALUE"] for _, p in b
    }


@pytest.mark.parametrize("missing", [0, 4])
def test_missing_outer_or_inner_reference_keeps_original(missing):
    original = landscape()
    original.pop(missing)
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert derived == original
    assert audit[0][1]["decision"] == "retain_missing_reference"


@pytest.mark.parametrize("status", ["mixed_spacing", "insufficient_levels"])
def test_uncertain_sheet_interval_prevents_recovery(status):
    original = landscape()
    derived, audit = recover_contour_heights(
        original, SHEET, {"status": status, "interval_m": None}
    )
    assert derived == original
    assert audit[0][1]["decision"] == "retain_uncertain_interval"


def test_inconsistent_inner_heights_or_same_outer_and_inner_are_ambiguous():
    original = landscape()
    original.extend(
        records([([(30, 30), (40, 30), (40, 40), (30, 40), (30, 30)], 160)])
    )
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert derived == original
    assert audit[0][1]["decision"] == "retain_ambiguous_nesting"
    original = landscape()
    original[4][1]["VALUE"] = 100
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert derived == original
    assert audit[0][1]["decision"] == "retain_ambiguous_nesting"


def test_neither_recorded_height_matches_the_bracket():
    original = landscape()
    original[1][1]["VALUE"] = 160
    original[2][1]["VALUE"] = original[3][1]["VALUE"] = 180
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert derived == original
    assert audit[0][1]["decision"] == "retain_neither_height_supported"


def test_candidate_near_boundary_stays_unresolved():
    original = landscape()
    derived, audit = recover_contour_heights(
        original, box(-200, -65, 200, 200), inference()
    )
    assert derived == original
    assert audit[0][1]["decision"] == "retain_near_sheet_boundary"


@pytest.mark.parametrize(
    "line,decision",
    [
        ([(60, 0), (80, 0)], "retain_contact_or_crossing"),
        ([(20, 0), (35, 0)], "retain_reference_contact_or_crossing"),
        ([(30, 0), (40, 0)], "retain_open_inner_reference"),
        ([(70, 0), (80, 0)], "retain_open_outer_reference"),
    ],
)
def test_contacts_and_incomplete_inner_linework_prevent_recovery(line, decision):
    original = landscape() + records([(line, 130)])
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert derived == original
    assert audit[0][1]["decision"] == decision


def test_off_series_inner_loops_are_preserved_and_cannot_contradict_relief():
    off_series = [(30, 30), (40, 30), (40, 40), (30, 40), (30, 30)]
    original = landscape() + records([(off_series, 130)])
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert audit[0][1]["decision"] == "supported_recovery"
    assert derived[-1] == original[-1]
    original[-1][1]["VALUE"] = 110
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert derived == original
    assert audit[0][1]["decision"] == "retain_contradictory_off_series_reference"


def test_nested_unresolved_loops_cannot_support_each_other():
    original = landscape() + records(
        [
            (square(40)[:4], 120),
            ([(40, -40), (-40, -40)], 140),
        ]
    )
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert derived == original
    assert len(audit) == 2
    assert {p["decision"] for _, p in audit} == {"retain_unresolved_neighbour"}


def test_off_series_enclosing_contour_must_fit_the_same_relief_direction():
    original = landscape() + records([(square(80), 110)])
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert audit[0][1]["decision"] == "supported_recovery"
    assert derived[-1] == original[-1]
    original[-1][1]["VALUE"] = 130
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert derived == original
    assert audit[0][1]["decision"] == "retain_contradictory_off_series_reference"


def test_four_arm_contact_lens_is_not_an_attribute_repair():
    original = records(
        [
            ([(-100, -50), (0, 0), (5, 5), (10, 0), (100, -50)], 1080),
            ([(-100, 50), (0, 0), (5, -5), (10, 0), (100, 50)], 1100),
        ]
    )
    assert recover_contour_heights(original, SHEET, inference()) == (original, [])


def tiny_fragment_loop(gap=6, width=30, target=880, split=False):
    lines = [([(gap, 0), (width, 0), (width, width), (0, width), (0, 0)], target)]
    fragments = (
        [[(0, 0), (gap / 2, 0)], [(gap / 2, 0), (gap, 0)]]
        if split
        else [[(0, 0), (gap, 0)]]
    )
    return records(lines + [(xy, 1240) for xy in fragments])


@pytest.mark.parametrize(
    "gap,width,accepted",
    [(6, 30, True), (10, 30, True), (10.1, 30, False), (6, 15, True), (6.1, 15, False)],
)
def test_fragment_repair_obeys_both_absolute_and_relative_length_caps(
    gap, width, accepted
):
    original = tiny_fragment_loop(gap, width)
    snapshot = copy.deepcopy(original)
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert (audit[0][1]["decision"] == "supported_recovery") == accepted
    assert original == snapshot
    if accepted:
        assert all(p["VALUE"] == 880 for _, p in derived)
        assert audit[0][1]["method"] == "small_fragment_between_same_height_endpoints"
        assert audit[0][1]["fragment_recovery"]["length_m"] == pytest.approx(gap)
        assert audit[0][1]["fragment_recovery"]["fraction_of_loop"] == pytest.approx(
            gap / (4 * width)
        )
        assert derived[1][1]["height_recovery"]["original_elevation_m"] == 1240
    else:
        assert derived == original


def test_splitting_an_overlong_fragment_does_not_bypass_the_cap():
    original = tiny_fragment_loop(12, 40, split=True)
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert derived == original
    assert audit[0][1]["fragment_recovery"]["length_m"] == 12
    assert audit[0][1]["fragment_recovery"]["decision"] == "retain_size_limit"
    original = tiny_fragment_loop(6, 30, split=True)
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert len([p for _, p in derived if "height_recovery" in p]) == 2
    assert audit[0][1]["changed_length_m"] == 6


@pytest.mark.parametrize(
    "limits", [{"max_fragment_length_m": 0}, {"max_fragment_fraction": 0}]
)
def test_fragment_fallback_can_be_disabled(limits):
    original = tiny_fragment_loop()
    derived, audit = recover_contour_heights(original, SHEET, inference(), **limits)
    assert derived == original
    assert audit[0][1]["fragment_recovery"]["decision"] == "disabled"


def test_small_fragment_still_rejects_branches_gaps_edges_and_off_series_heights():
    branch = tiny_fragment_loop() + records([([(3, 0), (3, -20)], 900)])
    assert recover_contour_heights(branch, SHEET, inference())[0] == branch
    gap = tiny_fragment_loop()
    gap[1] = (LineString([(-0.01, 0), (6, 0)]), gap[1][1])
    assert recover_contour_heights(gap, SHEET, inference()) == (gap, [])
    original = tiny_fragment_loop()
    assert (
        recover_contour_heights(original, box(-10, -10, 100, 100), inference())[0]
        == original
    )
    off_series = tiny_fragment_loop(target=885)
    assert recover_contour_heights(off_series, SHEET, inference())[0] == off_series


def test_small_fragment_does_not_override_a_contradictory_available_reference():
    original = tiny_fragment_loop() + records([(square(100), 920)])
    derived, audit = recover_contour_heights(original, SHEET, inference())
    assert derived == original
    assert (
        audit[0][1]["fragment_recovery"]["decision"] == "retain_contradictory_reference"
    )


@pytest.mark.parametrize(
    "limits",
    [
        {"max_fragment_length_m": -1},
        {"max_fragment_length_m": float("nan")},
        {"max_fragment_fraction": 0.5},
        {"max_fragment_fraction": float("inf")},
    ],
)
def test_invalid_fragment_limits_are_rejected(limits):
    with pytest.raises(ValueError, match="Height repair"):
        recover_contour_heights(tiny_fragment_loop(), SHEET, inference(), **limits)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
