"""Test the periodic function."""

import pytest

from microgen import Phase, Rve, periodic_split_and_translate, shape
from microgen.cad import CadShape

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/
# ruff: noqa: E501 line-too-long https://docs.astral.sh/ruff/rules/line-too-long/

N_PARTS_NO_INTERSECTION = 1
N_PARTS_ON_FACE = 2
N_PARTS_ON_EDGE = 4
N_PARTS_ON_CORNER = 8


def _generate_sphere(x: float, y: float, z: float, rve: Rve) -> Phase:
    """Generate a sphere at the given position and test periodicity."""
    elem = shape.sphere.Sphere(center=(x, y, z), radius=0.1)
    phase = Phase(shape=elem.generate())
    return periodic_split_and_translate(phase=phase, rve=rve)


def test_periodic_generates_warning_on_intersection_with_opposite_faces() -> None:
    """Test that x- and x+ faces of the RVE are intersected by the capsule."""
    rve = Rve(dim=1, center=(0.5, 0.5, 0.5))

    elem = shape.capsule.Capsule(center=(0.5, 0, 0.5), height=1, radius=0.1)
    phase = Phase(shape=elem.generate())

    expected_warning_msg = r"Object intersecting ([xyz])\+ and ([xyz])\- faces: not doing anything in this direction"
    with pytest.warns(UserWarning, match=expected_warning_msg):
        phase = periodic_split_and_translate(phase=phase, rve=rve)


def test_periodic_when_no_intersection() -> None:
    """Test that the periodic function does not raise an error when there is no intersection."""
    rve = Rve(dim=1, center=(0.5, 0.5, 0.5))
    phase = _generate_sphere(x=0.5, y=0.5, z=0.5, rve=rve)
    assert len(phase.solids) == N_PARTS_NO_INTERSECTION


@pytest.mark.parametrize(
    ("x", "y", "z"),
    [
        (0, 0.5, 0.5),
        (1, 0.5, 0.5),
        (0.5, 0, 0.5),
        (0.5, 1, 0.5),
        (0.5, 0.5, 0),
        (0.5, 0.5, 1),
    ],
)
def test_periodic_when_intersection_with_one_face(x: float, y: float, z: float) -> None:
    """Test that the periodic function does not raise an error when there is an intersection with one face."""
    rve = Rve(dim=1, center=(0.5, 0.5, 0.5))
    phase = _generate_sphere(x=x, y=y, z=z, rve=rve)
    assert len(phase.solids) == N_PARTS_ON_FACE


@pytest.mark.parametrize(
    ("x", "y", "z"),
    [
        (0, 0, 0.5),
        (0, 1, 0.5),
        (0, 0.5, 0),
        (0, 0.5, 1),
        (1, 0, 0.5),
        (1, 1, 0.5),
        (1, 0.5, 0),
        (1, 0.5, 1),
        (0.5, 0, 0),
        (0.5, 0, 1),
        (0.5, 1, 0),
        (0.5, 1, 1),
    ],
)
def test_periodic_when_intersection_with_one_edge(x: float, y: float, z: float) -> None:
    """Test that the periodic function does not raise an error when there is an intersection with one face."""
    rve = Rve(dim=1, center=(0.5, 0.5, 0.5))
    phase = _generate_sphere(x=x, y=y, z=z, rve=rve)
    assert len(phase.solids) == N_PARTS_ON_EDGE


@pytest.mark.parametrize(
    ("x", "y", "z"),
    [
        (0, 0, 0),
        (0, 0, 1),
        (0, 1, 0),
        (0, 1, 1),
        (1, 0, 0),
        (1, 0, 1),
        (1, 1, 0),
        (1, 1, 1),
    ],
)
def test_periodic_when_intersection_with_one_corner(
    x: float,
    y: float,
    z: float,
) -> None:
    """Test the periodic function."""
    rve = Rve(dim=1, center=(0.5, 0.5, 0.5))

    phase = _generate_sphere(x=x, y=y, z=z, rve=rve)
    assert len(phase.solids) == N_PARTS_ON_CORNER


# Volume- and bounds-based regression checks for the OCP-direct rewrite of
# microgen.periodic.  The split/translate path was rewritten from
# cadquery's ``Workplane.split() + .solids(">X")`` to
# ``BRepAlgoAPI_Splitter`` + a centroid-side selector keyed on
# ``_SIDE_DIR``.  A sign error there would still produce the right *count*
# of solids (the existing parametric tests above only check counts), so we
# pin two harder invariants:
#
#   1. every split fragment lands strictly inside the RVE box, and
#   2. the total fragment volume equals the original sphere volume
#      (mass conservation under the periodic rearrangement).
_VOL_REL_TOL = 5e-3  # OCCT booleans have small volumetric drift


def _total_volume(phase: Phase) -> float:
    # phase.solids is a list of raw TopoDS_Solid; wrap each so we can call
    # the CadShape helpers (Volume / BoundingBox).
    return sum(float(CadShape(s).Volume()) for s in phase.solids)


def _all_solids_inside(phase: Phase, rve: Rve, tol: float = 1e-3) -> bool:
    # OCCT boolean intersection (BRepAlgoAPI_Common with the RVE box) leaves a
    # small geometric drift on the order of 1e-4 in the resulting bounding
    # box.  A misplaced fragment from a sign error would be off by O(rve.dim),
    # so 1e-3 is loose enough for OCCT noise yet tight enough to catch the
    # bug class this test targets.
    for solid in phase.solids:
        bb = CadShape(solid).BoundingBox()
        if (
            bb.xmin < rve.min_point[0] - tol
            or bb.xmax > rve.max_point[0] + tol
            or bb.ymin < rve.min_point[1] - tol
            or bb.ymax > rve.max_point[1] + tol
            or bb.zmin < rve.min_point[2] - tol
            or bb.zmax > rve.max_point[2] + tol
        ):
            return False
    return True


@pytest.mark.parametrize(
    ("x", "y", "z", "expected_parts"),
    [
        # face crossings — one cut plane, two fragments
        (0.0, 0.5, 0.5, N_PARTS_ON_FACE),
        (0.5, 1.0, 0.5, N_PARTS_ON_FACE),
        # edge crossings — two cut planes, four fragments
        (0.0, 0.0, 0.5, N_PARTS_ON_EDGE),
        (1.0, 0.5, 1.0, N_PARTS_ON_EDGE),
        # corner crossings — three cut planes, eight fragments
        (0.0, 0.0, 0.0, N_PARTS_ON_CORNER),
        (1.0, 1.0, 1.0, N_PARTS_ON_CORNER),
    ],
)
def test_periodic_split_conserves_volume_and_stays_inside_rve(
    x: float,
    y: float,
    z: float,
    expected_parts: int,
) -> None:
    """Periodic split must conserve volume and land every fragment in the RVE.

    Catches sign / direction bugs in the rewritten ``_SIDE_DIR`` translate
    logic that the count-only parametric tests above cannot see.
    """
    import math

    rve = Rve(dim=1, center=(0.5, 0.5, 0.5))
    radius = 0.1
    expected_volume = (4.0 / 3.0) * math.pi * radius**3

    phase = _generate_sphere(x=x, y=y, z=z, rve=rve)

    assert len(phase.solids) == expected_parts
    assert _all_solids_inside(phase, rve), (
        f"At least one fragment lies outside the RVE for seed ({x}, {y}, {z}); "
        "this typically means a translate(...) call moved the wrong half."
    )
    assert _total_volume(phase) == pytest.approx(
        expected_volume,
        rel=_VOL_REL_TOL,
    ), (
        f"Total fragment volume drifted from the sphere volume for seed "
        f"({x}, {y}, {z}); this typically means a fragment was duplicated "
        f"or dropped by the split-and-translate path."
    )
