"""Tests for the F-rep implicit operations and Shape implicit capabilities."""

from __future__ import annotations

import numpy as np
import pytest
import pyvista as pv

from microgen.shape.implicit_ops import (
    batch_smooth_union,
    blend,
    difference,
    from_field,
    intersection,
    repeat,
    shell,
    smooth_difference,
    smooth_intersection,
    smooth_union,
    union,
)
from microgen.shape.shape import Shape


# ---------------------------------------------------------------------------
# Helpers — inline SDF lambdas (no primitive factories needed)
# ---------------------------------------------------------------------------


def _sphere_field(cx=0.0, cy=0.0, cz=0.0, r=1.0):
    """Return a sphere SDF function and bounds."""
    margin = r * 1.1
    return (
        lambda x, y, z: np.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2) - r,
        (cx - margin, cx + margin, cy - margin, cy + margin, cz - margin, cz + margin),
    )


def _box_field(cx=0.0, cy=0.0, cz=0.0, hx=0.5, hy=0.5, hz=0.5):
    """Return a box SDF function and bounds."""
    margin = max(hx, hy, hz) * 0.1

    def sdf(x, y, z):
        qx = np.abs(x - cx) - hx
        qy = np.abs(y - cy) - hy
        qz = np.abs(z - cz) - hz
        outside = np.sqrt(
            np.maximum(qx, 0.0) ** 2
            + np.maximum(qy, 0.0) ** 2
            + np.maximum(qz, 0.0) ** 2
        )
        inside = np.minimum(np.maximum(qx, np.maximum(qy, qz)), 0.0)
        return outside + inside

    return (
        sdf,
        (
            cx - hx - margin,
            cx + hx + margin,
            cy - hy - margin,
            cy + hy + margin,
            cz - hz - margin,
            cz + hz + margin,
        ),
    )


def _make_sphere(cx=0.0, cy=0.0, cz=0.0, r=1.0):
    func, bounds = _sphere_field(cx, cy, cz, r)
    return Shape(func=func, bounds=bounds)


def _make_box(cx=0.0, cy=0.0, cz=0.0, hx=0.5, hy=0.5, hz=0.5):
    func, bounds = _box_field(cx, cy, cz, hx, hy, hz)
    return Shape(func=func, bounds=bounds)


# ---------------------------------------------------------------------------
# Evaluate
# ---------------------------------------------------------------------------


class TestEvaluate:
    """Test implicit field evaluation."""

    def test_sphere_inside(self):
        s = _make_sphere()
        assert s.evaluate(np.array([0.0]), np.array([0.0]), np.array([0.0]))[0] < 0

    def test_sphere_surface(self):
        s = _make_sphere()
        val = s.evaluate(np.array([1.0]), np.array([0.0]), np.array([0.0]))[0]
        assert abs(val) < 1e-10

    def test_sphere_outside(self):
        s = _make_sphere()
        assert s.evaluate(np.array([2.0]), np.array([0.0]), np.array([0.0]))[0] > 0

    def test_no_func_raises(self):
        s = Shape()
        with pytest.raises(ValueError, match="No implicit scalar field"):
            s.evaluate(np.array([0.0]), np.array([0.0]), np.array([0.0]))


# ---------------------------------------------------------------------------
# Boolean operations
# ---------------------------------------------------------------------------


class TestBooleans:
    """Test hard boolean operations."""

    def test_union_min(self):
        s1 = _make_sphere(cx=-0.5)
        s2 = _make_sphere(cx=0.5)
        u = union(s1, s2)
        x, y, z = np.array([0.0]), np.array([0.0]), np.array([0.0])
        expected = min(s1.evaluate(x, y, z)[0], s2.evaluate(x, y, z)[0])
        assert u.evaluate(x, y, z)[0] == pytest.approx(expected)

    def test_intersection_max(self):
        s1 = _make_sphere(cx=-0.5)
        s2 = _make_sphere(cx=0.5)
        i = intersection(s1, s2)
        x, y, z = np.array([0.0]), np.array([0.0]), np.array([0.0])
        expected = max(s1.evaluate(x, y, z)[0], s2.evaluate(x, y, z)[0])
        assert i.evaluate(x, y, z)[0] == pytest.approx(expected)

    def test_difference(self):
        s1 = _make_sphere()
        s2 = _make_sphere(cx=0.5, r=0.5)
        d = difference(s1, s2)
        x, y, z = np.array([-0.5]), np.array([0.0]), np.array([0.0])
        assert d.evaluate(x, y, z)[0] < 0

    def test_operators(self):
        s1 = _make_sphere()
        s2 = _make_box()
        x, y, z = np.array([0.0]), np.array([0.0]), np.array([0.0])

        u = s1 | s2
        assert u.evaluate(x, y, z)[0] < 0

        i = s1 & s2
        assert i.evaluate(x, y, z)[0] < 0

        d = s1 - s2
        val = d.evaluate(x, y, z)[0]
        assert isinstance(val, (float, np.floating))

    def test_complement(self):
        s = _make_sphere()
        c = ~s
        x, y, z = np.array([0.0]), np.array([0.0]), np.array([0.0])
        assert c.evaluate(x, y, z)[0] > 0


# ---------------------------------------------------------------------------
# Smooth booleans
# ---------------------------------------------------------------------------


class TestSmoothBooleans:
    """Test smooth boolean operations."""

    def test_smooth_union_blending_zone(self):
        s1 = _make_sphere(cx=-0.5)
        s2 = _make_sphere(cx=0.5)
        su = smooth_union(s1, s2, k=0.5)
        hard = union(s1, s2)
        x, y, z = np.array([0.0]), np.array([0.0]), np.array([0.0])
        assert su.evaluate(x, y, z)[0] <= hard.evaluate(x, y, z)[0] + 1e-10

    def test_smooth_union_k_zero_equals_hard(self):
        s1 = _make_sphere(cx=-0.5)
        s2 = _make_sphere(cx=0.5)
        su = smooth_union(s1, s2, k=0.0)
        hard = union(s1, s2)
        x = np.linspace(-2, 2, 10)
        y = np.zeros(10)
        z = np.zeros(10)
        np.testing.assert_allclose(
            su.evaluate(x, y, z), hard.evaluate(x, y, z), atol=1e-10
        )

    def test_smooth_intersection(self):
        s1 = _make_sphere()
        s2 = _make_box()
        si = smooth_intersection(s1, s2, k=0.3)
        x, y, z = np.array([0.0]), np.array([0.0]), np.array([0.0])
        assert si.evaluate(x, y, z)[0] < 0

    def test_smooth_difference(self):
        s1 = _make_sphere()
        s2 = _make_sphere(cx=1.0, r=0.5)
        sd = smooth_difference(s1, s2, k=0.2)
        x, y, z = np.array([-0.5]), np.array([0.0]), np.array([0.0])
        assert sd.evaluate(x, y, z)[0] < 0

    def test_smooth_methods(self):
        s1 = _make_sphere()
        s2 = _make_box()
        x, y, z = np.array([0.0]), np.array([0.0]), np.array([0.0])

        su = s1.smooth_union(s2, k=0.3)
        assert isinstance(su.evaluate(x, y, z)[0], (float, np.floating))

        si = s1.smooth_intersection(s2, k=0.3)
        assert isinstance(si.evaluate(x, y, z)[0], (float, np.floating))

        sd = s1.smooth_difference(s2, k=0.3)
        assert isinstance(sd.evaluate(x, y, z)[0], (float, np.floating))


# ---------------------------------------------------------------------------
# Transforms
# ---------------------------------------------------------------------------


class TestTransforms:
    """Test implicit field transform operations."""

    def test_translate(self):
        s = _make_sphere()
        st = s.translate((2, 0, 0))
        x, y, z = np.array([2.0]), np.array([0.0]), np.array([0.0])
        assert st.evaluate(x, y, z)[0] < 0

        x0 = np.array([0.0])
        assert st.evaluate(x0, y, z)[0] > 0

    def test_rotate_90(self):
        # Elongated box along x, rotate 90 around z -> elongated along y
        func, bounds = _box_field(hx=1.0, hy=0.1, hz=0.1)
        box = Shape(func=func, bounds=bounds)
        rotated = box.rotate((0, 0, 90), convention="xyz")
        # Point along y axis should be inside
        assert (
            rotated.evaluate(np.array([0.0]), np.array([0.5]), np.array([0.0]))[0] < 0
        )
        # Point along x axis (was inside, now outside)
        assert (
            rotated.evaluate(np.array([0.5]), np.array([0.0]), np.array([0.0]))[0] > 0
        )

    def test_scale(self):
        s = _make_sphere()
        ss = s.scale(2.0)
        x, y, z = np.array([1.5]), np.array([0.0]), np.array([0.0])
        assert s.evaluate(x, y, z)[0] > 0
        assert ss.evaluate(x, y, z)[0] < 0

    def test_translate_bounds(self):
        s = _make_sphere()
        st = s.translate((5, 0, 0))
        assert st.bounds is not None
        assert st.bounds[0] > 3.0

    def test_scale_bounds(self):
        s = _make_sphere()
        ss = s.scale(3.0)
        assert ss.bounds is not None
        assert ss.bounds[1] > 3.0


# ---------------------------------------------------------------------------
# generate_surface_mesh
# ---------------------------------------------------------------------------


class TestGenerateVtk:
    """Test default mesh generation from implicit field."""

    def test_sphere_mesh(self):
        s = _make_sphere()
        mesh = s.generate_surface_mesh(resolution=30)
        assert isinstance(mesh, pv.PolyData)
        assert mesh.n_cells > 0

    def test_box_mesh(self):
        b = _make_box()
        mesh = b.generate_surface_mesh(resolution=30)
        assert isinstance(mesh, pv.PolyData)
        assert mesh.n_cells > 0

    def test_boolean_mesh(self):
        s1 = _make_sphere()
        s2 = _make_box(hx=0.6, hy=0.6, hz=0.6)
        result = s1 & s2
        mesh = result.generate_surface_mesh(resolution=30)
        assert isinstance(mesh, pv.PolyData)
        assert mesh.n_cells > 0

    def test_no_func_raises(self):
        s = Shape()
        with pytest.raises(NotImplementedError, match="No implicit field"):
            s.generate_surface_mesh()

    def test_no_bounds_raises(self):
        s = Shape(func=lambda x, y, z: x**2 + y**2 + z**2 - 1)
        with pytest.raises(ValueError, match="Bounds must be provided"):
            s.generate_surface_mesh()

    def test_explicit_bounds_override(self):
        s = _make_sphere()
        mesh = s.generate_surface_mesh(bounds=(-2, 2, -2, 2, -2, 2), resolution=30)
        assert mesh.n_cells > 0


# ---------------------------------------------------------------------------
# Bounds propagation
# ---------------------------------------------------------------------------


class TestBoundsPropagation:
    """Test bounds merge behavior."""

    def test_union_expands(self):
        s1 = _make_sphere(cx=-2)
        s2 = _make_sphere(cx=2)
        u = s1 | s2
        assert u.bounds[0] <= s1.bounds[0]
        assert u.bounds[1] >= s2.bounds[1]

    def test_intersection_shrinks(self):
        s1 = _make_sphere(cx=-0.5, r=1.5)
        s2 = _make_sphere(cx=0.5, r=1.5)
        i = s1 & s2
        assert i.bounds[0] >= s1.bounds[0]
        assert i.bounds[1] <= s2.bounds[1]

    def test_none_bounds(self):
        f = lambda x, y, z: x  # noqa: E731
        p = Shape(func=f, bounds=None)
        s = _make_sphere()
        u = p | s
        assert u.bounds == s.bounds


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


class TestUtilities:
    """Test utility functions."""

    def test_shell(self):
        s = _make_sphere()
        sh = shell(s, thickness=0.2)
        # At radius=1.0 (surface), |f|=0, shell should be inside
        assert sh.evaluate(np.array([1.0]), np.array([0.0]), np.array([0.0]))[0] < 0
        # At origin, |f|=1.0 >> 0.1, shell should be outside
        assert sh.evaluate(np.array([0.0]), np.array([0.0]), np.array([0.0]))[0] > 0

    def test_repeat(self):
        s = _make_sphere(r=0.3)
        r = repeat(s, spacing=(1.0, 1.0, 1.0))
        # At (1,0,0) should be inside a repeated copy
        assert r.evaluate(np.array([1.0]), np.array([0.0]), np.array([0.0]))[0] < 0
        assert r.bounds is None

    def test_blend(self):
        s1 = _make_sphere()
        s2 = _make_box()
        b = blend(s1, s2, factor=0.5)
        x, y, z = np.array([0.0]), np.array([0.0]), np.array([0.0])
        expected = 0.5 * s1.evaluate(x, y, z)[0] + 0.5 * s2.evaluate(x, y, z)[0]
        assert b.evaluate(x, y, z)[0] == pytest.approx(expected)

    def test_from_field(self):
        shape = from_field(
            func=lambda x, y, z: x**2 + y**2 + z**2 - 1,
            bounds=(-2, 2, -2, 2, -2, 2),
        )
        assert shape.evaluate(np.array([0.0]), np.array([0.0]), np.array([0.0]))[0] < 0
        assert shape.bounds == (-2, 2, -2, 2, -2, 2)

    def test_batch_smooth_union(self):
        spheres = [_make_sphere(cx=i * 0.5) for i in range(5)]
        combined = batch_smooth_union(spheres, k=0.1)
        x, y, z = np.array([0.0]), np.array([0.0]), np.array([0.0])
        assert combined.evaluate(x, y, z)[0] < 0

    def test_batch_smooth_union_empty_raises(self):
        with pytest.raises(ValueError, match="at least one shape"):
            batch_smooth_union([], k=0.1)

    def test_batch_smooth_union_hard(self):
        spheres = [_make_sphere(cx=i * 0.5) for i in range(3)]
        combined = batch_smooth_union(spheres, k=0.0)
        x, y, z = np.array([0.0]), np.array([0.0]), np.array([0.0])
        assert combined.evaluate(x, y, z)[0] < 0


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


class TestErrorHandling:
    """Test error conditions."""

    def test_transform_without_func_raises(self):
        s = Shape()
        with pytest.raises(ValueError, match="No implicit scalar field"):
            s.translate((1, 0, 0))
        with pytest.raises(ValueError, match="No implicit scalar field"):
            s.rotate((0, 0, 45))
        with pytest.raises(ValueError, match="No implicit scalar field"):
            s.scale(2.0)

    def test_boolean_without_func_raises(self):
        s = Shape()
        other = _make_sphere()
        with pytest.raises(ValueError, match="No implicit scalar field"):
            s | other  # noqa: B015
        with pytest.raises(ValueError, match="No implicit scalar field"):
            ~s  # noqa: B015
