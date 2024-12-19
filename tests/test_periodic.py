"""Test the periodic function."""

import pytest

from microgen import Phase, Rve, periodic_split_and_translate, shape

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
