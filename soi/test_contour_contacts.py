# /// script
# requires-python = ">=3.12"
# dependencies = ["shapely>=2.1", "pytest>=8"]
# ///
import pytest
from shapely.geometry import LineString, Point, box

from contour_contacts import relax_contour_contacts

SHEET = box(-1000, -1000, 1000, 1000)


def records(lines):
    return [
        (LineString(xy), {"VALUE": z, "source_id": str(i)})
        for i, (xy, z) in enumerate(lines)
    ]


def test_tangent_contours_keep_heights_and_all_geometry_outside_contact():
    original = records(
        [
            ([(-100, 100), (0, 0), (100, 100)], 80),
            ([(-100, -100), (0, 0), (100, -100)], 120),
        ]
    )
    kept, contacts, relaxed = relax_contour_contacts(original, 15, SHEET)
    assert len(contacts) == 1
    assert contacts[0][1]["geometry_pattern"] == "tangent_contact"
    assert {p["VALUE"] for _, p in kept} == {80, 120}
    assert all(line.distance(Point(0, 0)) >= 14.9 for line, _ in kept)
    assert sum(g.length for g, _ in kept + relaxed) == pytest.approx(
        sum(g.length for g, _ in original)
    )
    assert all(len(g.coords) == 3 for g, _ in original)


@pytest.mark.parametrize("split", [False, True])
def test_small_contact_lens_works_with_split_source_lines(split):
    lines = [
        ([(-100, -50), (0, 0), (5, 5), (10, 0), (100, -50)], 1080),
        ([(-100, 50), (0, 0), (5, -5), (10, 0), (100, 50)], 1100),
    ]
    if split:
        lines = [(part, z) for xy, z in lines for part in [xy[:2], xy[1:4], xy[3:]]]
    original = records(lines)
    kept, contacts, relaxed = relax_contour_contacts(original, 15, SHEET)
    assert len(contacts) == 2
    assert {p["geometry_pattern"] for _, p in contacts} == {"short_contact_lens"}
    assert {p["VALUE"] for _, p in kept} == {1080, 1100}
    assert sum(g.length for g, _ in kept + relaxed) == pytest.approx(
        sum(g.length for g, _ in original)
    )
    assert all(
        line.distance(point) >= 14.9 for line, _ in kept for point, _ in contacts
    )


@pytest.mark.parametrize(
    "lines",
    [
        [
            ([(-100, -100), (0, 0), (100, 100)], 80),
            ([(-100, 100), (0, 0), (100, -100)], 100),
        ],  # Isolated crossing.
        [
            ([(-100, 0), (0, 0), (100, 0)], 80),
            ([(0, 0), (0, 100)], 100),
        ],  # Form-like T-junction.
        [([(-100, 0), (0, 0)], 80), ([(0, 0), (100, 0)], 100)],  # Height change.
        [
            ([(-100, 0), (0, 0), (100, 0)], 80),
            ([(-100, 0), (0, 0), (100, 0)], 100),
        ],  # Shared segment.
        [
            ([(-100, -50), (0, 0), (5, 100), (10, 0), (100, -50)], 80),
            ([(-100, 50), (0, 0), (5, -100), (10, 0), (100, 50)], 100),
        ],  # Large loop.
    ],
)
def test_other_conflicts_are_not_relaxed(lines):
    original = records(lines)
    kept, contacts, relaxed = relax_contour_contacts(original, 15, SHEET)
    assert contacts == relaxed == []
    assert kept == original


def test_boundary_contacts_and_wholly_enclosed_contours_remain_strict():
    tangent = records(
        [
            ([(-100, 100), (0, 0), (100, 100)], 80),
            ([(-100, -100), (0, 0), (100, -100)], 100),
        ]
    )
    assert relax_contour_contacts(tangent, 15, box(-5, -1000, 1000, 1000))[1] == []
    tiny_ring = records(
        [
            ([(0, 0), (1, 1), (2, 0), (1, -1), (0, 0)], 100),
            ([(-100, -100), (0, 0), (-100, 100)], 80),
        ]
    )
    assert relax_contour_contacts(tiny_ring, 15, SHEET) == (tiny_ring, [], [])
    assert relax_contour_contacts(tangent, 0, SHEET) == (tangent, [], [])
    with pytest.raises(ValueError, match="radius"):
        relax_contour_contacts(tangent, -1, SHEET)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
