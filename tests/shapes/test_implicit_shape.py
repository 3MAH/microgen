"""Tests for the F-rep implicit modeling module."""

from __future__ import annotations

import numpy as np
import pytest
import pyvista as pv

from microgen.shape import implicit_shape, surface_functions, Tpms
from microgen.shape.implicit_shape import (
    ImplicitShape,
    blend,
    difference,
    from_field,
    implicit_box,
    implicit_capsule,
    implicit_cylinder,
    implicit_plane,
    implicit_sphere,
    implicit_torus,
    intersection,
    repeat,
    shell,
    smooth_difference,
    smooth_intersection,
    smooth_union,
    union,
)


# ---------------------------------------------------------------------------
# Primitive evaluation tests
# ---------------------------------------------------------------------------


class TestPrimitives:
    """Test that primitives evaluate correctly at known points."""

    def test_sphere_inside(self):
        s = implicit_sphere(center=(0, 0, 0), radius=1.0)
        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([0.0])
        assert s.evaluate(x, y, z)[0] < 0

    def test_sphere_surface(self):
        s = implicit_sphere(center=(0, 0, 0), radius=1.0)
        x = np.array([1.0])
        y = np.array([0.0])
        z = np.array([0.0])
        assert abs(s.evaluate(x, y, z)[0]) < 1e-10

    def test_sphere_outside(self):
        s = implicit_sphere(center=(0, 0, 0), radius=1.0)
        x = np.array([2.0])
        y = np.array([0.0])
        z = np.array([0.0])
        assert s.evaluate(x, y, z)[0] > 0

    def test_box_inside(self):
        b = implicit_box(center=(0, 0, 0), half_extents=(1, 1, 1))
        assert b.evaluate(np.array([0.0]), np.array([0.0]), np.array([0.0]))[0] < 0

    def test_box_outside(self):
        b = implicit_box(center=(0, 0, 0), half_extents=(1, 1, 1))
        assert b.evaluate(np.array([2.0]), np.array([0.0]), np.array([0.0]))[0] > 0

    def test_cylinder_inside(self):
        c = implicit_cylinder(center=(0, 0, 0), axis=(0, 0, 1), radius=1, height=2)
        assert c.evaluate(np.array([0.0]), np.array([0.0]), np.array([0.0]))[0] < 0

    def test_cylinder_outside(self):
        c = implicit_cylinder(center=(0, 0, 0), axis=(0, 0, 1), radius=1, height=2)
        assert c.evaluate(np.array([2.0]), np.array([0.0]), np.array([0.0]))[0] > 0

    def test_plane(self):
        p = implicit_plane(point=(0, 0, 0), normal=(0, 0, 1))
        assert p.evaluate(np.array([0.0]), np.array([0.0]), np.array([-1.0]))[0] < 0
        assert p.evaluate(np.array([0.0]), np.array([0.0]), np.array([1.0]))[0] > 0

    def test_plane_bounds_none(self):
        p = implicit_plane()
        assert p._bounds is None

    def test_torus_inside(self):
        t = implicit_torus(center=(0, 0, 0), major_r=1.0, minor_r=0.25)
        # Point on the ring center at (1, 0, 0) should be inside
        assert t.evaluate(np.array([1.0]), np.array([0.0]), np.array([0.0]))[0] < 0

    def test_torus_outside(self):
        t = implicit_torus(center=(0, 0, 0), major_r=1.0, minor_r=0.25)
        assert t.evaluate(np.array([0.0]), np.array([0.0]), np.array([0.0]))[0] > 0

    def test_capsule_inside(self):
        c = implicit_capsule(start=(0, 0, -1), end=(0, 0, 1), radius=0.5)
        assert c.evaluate(np.array([0.0]), np.array([0.0]), np.array([0.0]))[0] < 0

    def test_capsule_outside(self):
        c = implicit_capsule(start=(0, 0, -1), end=(0, 0, 1), radius=0.5)
        assert c.evaluate(np.array([2.0]), np.array([0.0]), np.array([0.0]))[0] > 0


# ---------------------------------------------------------------------------
# Boolean operations
# ---------------------------------------------------------------------------


class TestBooleans:
    """Test hard boolean operations."""

    def test_union_min(self):
        s1 = implicit_sphere(center=(-0.5, 0, 0), radius=1.0)
        s2 = implicit_sphere(center=(0.5, 0, 0), radius=1.0)
        u = union(s1, s2)
        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([0.0])
        expected = min(s1.evaluate(x, y, z)[0], s2.evaluate(x, y, z)[0])
        assert u.evaluate(x, y, z)[0] == pytest.approx(expected)

    def test_intersection_max(self):
        s1 = implicit_sphere(center=(-0.5, 0, 0), radius=1.0)
        s2 = implicit_sphere(center=(0.5, 0, 0), radius=1.0)
        i = intersection(s1, s2)
        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([0.0])
        expected = max(s1.evaluate(x, y, z)[0], s2.evaluate(x, y, z)[0])
        assert i.evaluate(x, y, z)[0] == pytest.approx(expected)

    def test_difference(self):
        s1 = implicit_sphere(center=(0, 0, 0), radius=1.0)
        s2 = implicit_sphere(center=(0.5, 0, 0), radius=0.5)
        d = difference(s1, s2)
        # Inside s1 but outside s2
        x = np.array([-0.5])
        y = np.array([0.0])
        z = np.array([0.0])
        assert d.evaluate(x, y, z)[0] < 0

    def test_operators(self):
        s1 = implicit_sphere(radius=1.0)
        s2 = implicit_box(half_extents=(0.5, 0.5, 0.5))
        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([0.0])

        u = s1 | s2
        assert u.evaluate(x, y, z)[0] < 0

        i = s1 & s2
        assert i.evaluate(x, y, z)[0] < 0

        d = s1 - s2
        # At origin, both are inside, so difference should be max(f_s1, -f_s2)
        val = d.evaluate(x, y, z)[0]
        assert isinstance(val, (float, np.floating))

    def test_complement(self):
        s = implicit_sphere(radius=1.0)
        c = ~s
        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([0.0])
        assert c.evaluate(x, y, z)[0] > 0  # origin was inside, complement is outside


# ---------------------------------------------------------------------------
# Smooth booleans
# ---------------------------------------------------------------------------


class TestSmoothBooleans:
    """Test smooth boolean operations."""

    def test_smooth_union_blending_zone(self):
        s1 = implicit_sphere(center=(-0.5, 0, 0), radius=1.0)
        s2 = implicit_sphere(center=(0.5, 0, 0), radius=1.0)
        su = smooth_union(s1, s2, k=0.5)
        hard = union(s1, s2)
        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([0.0])
        # Smooth union should be <= hard union (more material in blend zone)
        assert su.evaluate(x, y, z)[0] <= hard.evaluate(x, y, z)[0] + 1e-10

    def test_smooth_union_k_zero_equals_hard(self):
        s1 = implicit_sphere(center=(-0.5, 0, 0), radius=1.0)
        s2 = implicit_sphere(center=(0.5, 0, 0), radius=1.0)
        su = smooth_union(s1, s2, k=0.0)
        hard = union(s1, s2)
        x = np.linspace(-2, 2, 10)
        y = np.zeros(10)
        z = np.zeros(10)
        np.testing.assert_allclose(
            su.evaluate(x, y, z), hard.evaluate(x, y, z), atol=1e-10
        )

    def test_smooth_intersection(self):
        s1 = implicit_sphere(radius=1.0)
        s2 = implicit_box(half_extents=(0.5, 0.5, 0.5))
        si = smooth_intersection(s1, s2, k=0.3)
        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([0.0])
        assert si.evaluate(x, y, z)[0] < 0  # still inside

    def test_smooth_difference(self):
        s1 = implicit_sphere(radius=1.0)
        s2 = implicit_sphere(center=(1.0, 0, 0), radius=0.5)
        sd = smooth_difference(s1, s2, k=0.2)
        x = np.array([-0.5])
        y = np.array([0.0])
        z = np.array([0.0])
        assert sd.evaluate(x, y, z)[0] < 0  # far from subtraction

    def test_smooth_methods(self):
        s1 = implicit_sphere(radius=1.0)
        s2 = implicit_box(half_extents=(0.5, 0.5, 0.5))
        su = s1.smooth_union(s2, k=0.3)
        si = s1.smooth_intersection(s2, k=0.3)
        sd = s1.smooth_difference(s2, k=0.3)
        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([0.0])
        assert isinstance(su.evaluate(x, y, z)[0], (float, np.floating))
        assert isinstance(si.evaluate(x, y, z)[0], (float, np.floating))
        assert isinstance(sd.evaluate(x, y, z)[0], (float, np.floating))


# ---------------------------------------------------------------------------
# Transforms
# ---------------------------------------------------------------------------


class TestTransforms:
    """Test transform operations."""

    def test_translate(self):
        s = implicit_sphere(center=(0, 0, 0), radius=1.0)
        st = s.translate((2, 0, 0))
        x = np.array([2.0])
        y = np.array([0.0])
        z = np.array([0.0])
        assert st.evaluate(x, y, z)[0] < 0  # center moved to (2,0,0)

        x0 = np.array([0.0])
        assert st.evaluate(x0, y, z)[0] > 0  # origin is now outside

    def test_rotate_90(self):
        # Cylinder along z, rotate 90 degrees around y -> cylinder along x
        c = implicit_cylinder(radius=0.5, height=2.0)
        cr = c.rotate((0, 90, 0), convention="xyz")
        # Point along x axis (was z) should be inside
        assert cr.evaluate(np.array([0.5]), np.array([0.0]), np.array([0.0]))[0] < 0

    def test_scale(self):
        s = implicit_sphere(radius=1.0)
        ss = s.scale(2.0)
        # At distance 1.5, original sphere is outside but scaled should be inside
        x = np.array([1.5])
        y = np.array([0.0])
        z = np.array([0.0])
        assert s.evaluate(x, y, z)[0] > 0
        assert ss.evaluate(x, y, z)[0] < 0

    def test_translate_bounds(self):
        s = implicit_sphere(radius=1.0)
        st = s.translate((5, 0, 0))
        assert st._bounds is not None
        assert st._bounds[0] > 3.0  # xmin shifted

    def test_scale_bounds(self):
        s = implicit_sphere(radius=1.0)
        ss = s.scale(3.0)
        assert ss._bounds is not None
        assert ss._bounds[1] > 3.0  # xmax scaled


# ---------------------------------------------------------------------------
# generate_vtk
# ---------------------------------------------------------------------------


class TestGenerateVtk:
    """Test mesh generation."""

    def test_sphere_mesh(self):
        s = implicit_sphere(radius=1.0)
        mesh = s.generate_vtk(resolution=30)
        assert isinstance(mesh, pv.PolyData)
        assert mesh.n_cells > 0

    def test_box_mesh(self):
        b = implicit_box(half_extents=(0.5, 0.5, 0.5))
        mesh = b.generate_vtk(resolution=30)
        assert isinstance(mesh, pv.PolyData)
        assert mesh.n_cells > 0

    def test_boolean_mesh(self):
        s1 = implicit_sphere(radius=1.0)
        s2 = implicit_box(half_extents=(0.6, 0.6, 0.6))
        result = s1 & s2
        mesh = result.generate_vtk(resolution=30)
        assert isinstance(mesh, pv.PolyData)
        assert mesh.n_cells > 0

    def test_no_bounds_raises(self):
        shape = ImplicitShape(func=lambda x, y, z: x**2 + y**2 + z**2 - 1)
        with pytest.raises(ValueError, match="Bounds must be provided"):
            shape.generate_vtk()

    def test_explicit_bounds_override(self):
        s = implicit_sphere(radius=1.0)
        mesh = s.generate_vtk(bounds=(-2, 2, -2, 2, -2, 2), resolution=30)
        assert mesh.n_cells > 0

    def test_generateVtk_deprecated(self):
        s = implicit_sphere(radius=1.0)
        mesh = s.generateVtk(resolution=20)
        assert isinstance(mesh, pv.PolyData)


# ---------------------------------------------------------------------------
# Bounds propagation
# ---------------------------------------------------------------------------


class TestBoundsPropagation:
    """Test bounds merge behavior."""

    def test_union_expands(self):
        s1 = implicit_sphere(center=(-2, 0, 0), radius=1.0)
        s2 = implicit_sphere(center=(2, 0, 0), radius=1.0)
        u = s1 | s2
        assert u._bounds[0] <= s1._bounds[0]  # xmin from s1
        assert u._bounds[0] <= s2._bounds[0]
        assert u._bounds[1] >= s1._bounds[1]  # xmax includes both
        assert u._bounds[1] >= s2._bounds[1]

    def test_intersection_shrinks(self):
        s1 = implicit_sphere(center=(-0.5, 0, 0), radius=1.5)
        s2 = implicit_sphere(center=(0.5, 0, 0), radius=1.5)
        i = s1 & s2
        assert i._bounds[0] >= s1._bounds[0]
        assert i._bounds[1] <= s2._bounds[1]

    def test_plane_none_bounds(self):
        p = implicit_plane()
        s = implicit_sphere(radius=1.0)
        u = p | s
        # When one is None, result takes the other
        assert u._bounds == s._bounds


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


class TestUtilities:
    """Test utility functions."""

    def test_shell(self):
        s = implicit_sphere(radius=1.0)
        sh = shell(s, thickness=0.2)
        # At radius=1.0 (surface), |f|=0, shell should be inside (-thickness/2)
        x = np.array([1.0])
        y = np.array([0.0])
        z = np.array([0.0])
        assert sh.evaluate(x, y, z)[0] < 0
        # At origin, |f|=1.0 >> 0.1, shell should be outside
        x0 = np.array([0.0])
        assert sh.evaluate(x0, y, z)[0] > 0

    def test_repeat(self):
        s = implicit_sphere(radius=0.3)
        r = repeat(s, spacing=(1.0, 1.0, 1.0))
        # At (1,0,0) should be inside a repeated copy
        x = np.array([1.0])
        y = np.array([0.0])
        z = np.array([0.0])
        assert r.evaluate(x, y, z)[0] < 0
        assert r._bounds is None  # infinite repetition

    def test_blend(self):
        s1 = implicit_sphere(radius=1.0)
        s2 = implicit_box(half_extents=(0.5, 0.5, 0.5))
        b = blend(s1, s2, factor=0.5)
        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([0.0])
        expected = 0.5 * s1.evaluate(x, y, z)[0] + 0.5 * s2.evaluate(x, y, z)[0]
        assert b.evaluate(x, y, z)[0] == pytest.approx(expected)

    def test_from_field(self):
        shape = from_field(
            func=lambda x, y, z: x**2 + y**2 + z**2 - 1,
            bounds=(-2, 2, -2, 2, -2, 2),
        )
        assert shape.evaluate(np.array([0.0]), np.array([0.0]), np.array([0.0]))[0] < 0
        assert shape._bounds == (-2, 2, -2, 2, -2, 2)


# ---------------------------------------------------------------------------
# TPMS integration
# ---------------------------------------------------------------------------


class TestTpmsIntegration:
    """Test that TPMS still works and integrates with ImplicitShape."""

    def test_isinstance(self):
        tpms = Tpms(
            surface_function=surface_functions.gyroid,
            offset=0.3,
        )
        assert isinstance(tpms, ImplicitShape)

    def test_tpms_has_func(self):
        tpms = Tpms(
            surface_function=surface_functions.gyroid,
            offset=0.3,
        )
        assert tpms._func is not None

    def test_tpms_has_bounds(self):
        tpms = Tpms(
            surface_function=surface_functions.gyroid,
            offset=0.3,
        )
        assert tpms._bounds is not None

    def test_tpms_evaluate(self):
        tpms = Tpms(
            surface_function=surface_functions.gyroid,
            offset=0.3,
        )
        x = np.array([0.0])
        y = np.array([0.0])
        z = np.array([0.0])
        val = tpms.evaluate(x, y, z)
        assert isinstance(val[0], (float, np.floating))

    def test_tpms_boolean_with_sphere(self):
        tpms = Tpms(
            surface_function=surface_functions.gyroid,
            offset=0.3,
        )
        sphere = implicit_sphere(radius=0.4)
        result = tpms & sphere
        assert isinstance(result, ImplicitShape)
        mesh = result.generate_vtk(
            bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
            resolution=30,
        )
        assert isinstance(mesh, pv.PolyData)

    def test_tpms_generate_vtk_unchanged(self):
        """TPMS generate_vtk still works with its own override."""
        tpms = Tpms(
            surface_function=surface_functions.gyroid,
            offset=0.3,
        )
        mesh = tpms.generate_vtk(type_part="sheet")
        assert isinstance(mesh, pv.PolyData)
        assert mesh.n_cells > 0


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


class TestErrorHandling:
    """Test error conditions."""

    def test_no_func_evaluate_raises(self):
        shape = ImplicitShape()
        with pytest.raises(ValueError, match="No scalar field"):
            shape.evaluate(np.array([0.0]), np.array([0.0]), np.array([0.0]))

    def test_no_bounds_generate_vtk_raises(self):
        shape = ImplicitShape(func=lambda x, y, z: x)
        with pytest.raises(ValueError, match="Bounds must be provided"):
            shape.generate_vtk()
