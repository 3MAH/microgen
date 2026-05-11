"""Tests for the implicit (F-rep) field set on primitive shape classes.

Each primitive (`Box`, `Sphere`, `Cylinder`, `Capsule`, `Ellipsoid`) sets
``_func`` and ``_bounds`` in ``__init__``. These tests verify the SDF sign
inside / on / outside, AABB validity, and composability via the
operators.
"""

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/

from __future__ import annotations

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from microgen import Box, Capsule, Cylinder, Ellipsoid, Sphere


def _eval_scalar(shape, point):
    """Evaluate the implicit field at a single (x, y, z) point."""
    x = np.array([point[0]], dtype=np.float64)
    y = np.array([point[1]], dtype=np.float64)
    z = np.array([point[2]], dtype=np.float64)
    return float(shape.evaluate(x, y, z)[0])


@pytest.mark.parametrize(
    ("ctor", "kwargs", "inside_pt", "outside_pt", "surface_pt"),
    [
        # Sphere centered at origin, r=1
        (Sphere, dict(radius=1.0), (0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
        # Box (2, 1, 0.5) centered at origin
        (
            Box,
            dict(dim=(2.0, 1.0, 0.5)),
            (0.0, 0.0, 0.0),
            (2.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
        ),
        # Cylinder along x, r=0.5, h=2
        (
            Cylinder,
            dict(radius=0.5, height=2.0),
            (0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.5, 0.0),
        ),
        # Capsule along x, r=0.3, h=1
        (
            Capsule,
            dict(radius=0.3, height=1.0),
            (0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.3, 0.0),
        ),
        # Ellipsoid (1, 0.5, 0.25)
        (
            Ellipsoid,
            dict(radii=(1.0, 0.5, 0.25)),
            (0.0, 0.0, 0.0),
            (2.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
        ),
    ],
)
def test_field_sign(ctor, kwargs, inside_pt, outside_pt, surface_pt):
    """Inside negative, outside positive, surface ~= 0."""
    s = ctor(**kwargs)
    assert _eval_scalar(s, inside_pt) < 0.0
    assert _eval_scalar(s, outside_pt) > 0.0
    assert abs(_eval_scalar(s, surface_pt)) < 1e-9


@pytest.mark.parametrize(
    ("ctor", "kwargs"),
    [
        (Sphere, dict(radius=1.0)),
        (Box, dict(dim=(2.0, 1.0, 0.5))),
        (Cylinder, dict(radius=0.5, height=2.0)),
        (Capsule, dict(radius=0.3, height=1.0)),
        (Ellipsoid, dict(radii=(1.0, 0.5, 0.25))),
    ],
)
def test_aabb_contains_zero_set(ctor, kwargs):
    """Bounds must enclose the field's zero level set."""
    s = ctor(**kwargs)
    xmin, xmax, ymin, ymax, zmin, zmax = s.bounds
    # corners of the AABB are outside (field > 0); the center is inside (< 0)
    cx, cy, cz = (
        0.5 * (xmin + xmax),
        0.5 * (ymin + ymax),
        0.5 * (zmin + zmax),
    )
    assert _eval_scalar(s, (cx, cy, cz)) < 0.0
    for corner in [
        (xmin, ymin, zmin),
        (xmax, ymax, zmax),
        (xmin, ymax, zmin),
        (xmax, ymin, zmax),
    ]:
        assert _eval_scalar(s, corner) > 0.0


@pytest.mark.parametrize(
    ("ctor", "kwargs"),
    [
        (Sphere, dict(radius=1.0)),
        (Box, dict(dim=(2.0, 1.0, 0.5))),
        (Cylinder, dict(radius=0.5, height=2.0)),
        (Capsule, dict(radius=0.3, height=1.0)),
        (Ellipsoid, dict(radii=(1.0, 0.5, 0.25))),
    ],
)
def test_center_offset_translates_field(ctor, kwargs):
    """A non-zero ``center`` should translate the inside region accordingly."""
    offset = (3.0, -1.0, 2.0)
    s_origin = ctor(**kwargs)
    s_offset = ctor(**kwargs, center=offset)
    # The offset center is inside; the original origin is outside the offset shape.
    assert _eval_scalar(s_offset, offset) == pytest.approx(
        _eval_scalar(s_origin, (0.0, 0.0, 0.0)),
    )
    assert _eval_scalar(s_offset, (0.0, 0.0, 0.0)) > 0.0


def test_orientation_rotates_box_field():
    """Rotating a Box 90° about z swaps its x and y extents."""
    box = Box(dim=(2.0, 0.5, 0.5))
    rotated = Box(
        dim=(2.0, 0.5, 0.5),
        orientation=Rotation.from_euler("z", 90, degrees=True),
    )
    # Point along +y at distance 0.9 is outside the un-rotated box, inside the rotated one.
    assert _eval_scalar(box, (0.0, 0.9, 0.0)) > 0.0
    assert _eval_scalar(rotated, (0.0, 0.9, 0.0)) < 0.0


def test_union_field_is_min():
    """``a | b`` evaluates to ``min(f_a, f_b)`` pointwise."""
    a = Sphere(radius=0.5, center=(-0.6, 0.0, 0.0))
    b = Sphere(radius=0.5, center=(0.6, 0.0, 0.0))
    union = a | b
    # Each sphere's center is inside the union.
    assert _eval_scalar(union, (-0.6, 0.0, 0.0)) < 0.0
    assert _eval_scalar(union, (0.6, 0.0, 0.0)) < 0.0
    # The midpoint between them is outside both individual spheres
    # (radius 0.5, distance 0.6) and so outside the union.
    assert _eval_scalar(union, (0.0, 0.0, 0.0)) > 0.0


def test_intersection_field_is_max():
    """``a & b`` evaluates to ``max(f_a, f_b)`` pointwise."""
    sph = Sphere(radius=1.0)
    box = Box(dim=(1.0, 1.0, 1.0))
    inter = sph & box
    # Origin is inside both, so inside the intersection.
    assert _eval_scalar(inter, (0.0, 0.0, 0.0)) < 0.0
    # Point at (0.7, 0.0, 0.0) is inside the sphere but outside the box.
    assert _eval_scalar(sph, (0.7, 0.0, 0.0)) < 0.0
    assert _eval_scalar(box, (0.7, 0.0, 0.0)) > 0.0
    assert _eval_scalar(inter, (0.7, 0.0, 0.0)) > 0.0


def test_difference_field_carves_out_b():
    """``a - b`` removes b's interior from a."""
    sph = Sphere(radius=1.0)
    hole = Sphere(radius=0.3, center=(0.5, 0.0, 0.0))
    diff = sph - hole
    # Point at hole center: inside sph, inside hole → carved out.
    assert _eval_scalar(diff, (0.5, 0.0, 0.0)) > 0.0
    # Point at sphere center: inside sph, outside hole → kept.
    assert _eval_scalar(diff, (0.0, 0.0, 0.0)) < 0.0
