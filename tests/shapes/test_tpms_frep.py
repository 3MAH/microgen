"""Tests for TPMS F-rep integration with SDF normalization."""

from __future__ import annotations

import numpy as np
import pytest
import pyvista as pv

from microgen.shape import Tpms, surface_functions
from microgen.shape.implicit_ops import from_field, normalize_to_sdf, shell
from microgen.shape.shape import Shape


# ---------------------------------------------------------------------------
# SDF normalization quality
# ---------------------------------------------------------------------------


class TestSdfNormalization:
    """Test that SDF normalization produces well-behaved fields."""

    def test_gradient_magnitude_near_one(self):
        """After normalization, |grad(sdf)| should be approximately 1."""
        raw = from_field(surface_functions.gyroid)
        sdf_shape = normalize_to_sdf(raw)
        sdf = sdf_shape.field

        rng = np.random.default_rng(42)
        x = rng.uniform(-2, 2, 500)
        y = rng.uniform(-2, 2, 500)
        z = rng.uniform(-2, 2, 500)

        h = 1e-5
        gx = (sdf(x + h, y, z) - sdf(x - h, y, z)) / (2 * h)
        gy = (sdf(x, y + h, z) - sdf(x, y - h, z)) / (2 * h)
        gz = (sdf(x, y, z + h) - sdf(x, y, z - h)) / (2 * h)
        grad_mag = np.sqrt(gx**2 + gy**2 + gz**2)

        # Near the zero level set, gradient magnitude should be close to 1
        near_surface = np.abs(sdf(x, y, z)) < 0.3
        if near_surface.sum() > 10:
            assert np.mean(grad_mag[near_surface]) == pytest.approx(1.0, abs=0.3)

    def test_saddle_point_safety(self):
        """SDF at known saddle points should not produce NaN or Inf."""
        tpms = Tpms(
            surface_function=surface_functions.gyroid,
            offset=0.3,
        )
        # Gyroid saddle points are at (0, 0, 0) and equivalents
        x = np.array([0.0, 0.25, 0.5])
        y = np.array([0.0, 0.25, 0.5])
        z = np.array([0.0, 0.25, 0.5])
        vals = tpms.evaluate(x, y, z)
        assert np.all(np.isfinite(vals))

    def test_normalize_to_sdf_fallback(self):
        """normalize_to_sdf falls back to FD when autograd fails."""
        # A function using regular numpy (not autograd.numpy)
        import numpy as regular_np

        def non_autograd_field(x, y, z):
            return regular_np.sin(x) + regular_np.cos(y)

        shape = from_field(non_autograd_field, bounds=(-2, 2, -2, 2, -2, 2))
        sdf_shape = normalize_to_sdf(shape)
        vals = sdf_shape.field(np.array([0.0]), np.array([0.0]), np.array([0.0]))
        assert np.isfinite(vals[0])


# ---------------------------------------------------------------------------
# Tpms F-rep field
# ---------------------------------------------------------------------------


class TestTpmsFrepField:
    """Test that Tpms has a working F-rep implicit field."""

    def test_tpms_has_func(self):
        tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.3)
        assert tpms.field is not None

    def test_tpms_has_bounds(self):
        tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.3)
        assert tpms.bounds is not None

    def test_tpms_evaluate(self):
        tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.3)
        x = np.array([0.0, 0.1, 0.2])
        y = np.array([0.0, 0.1, 0.2])
        z = np.array([0.0, 0.1, 0.2])
        vals = tpms.evaluate(x, y, z)
        assert vals.shape == (3,)
        assert np.all(np.isfinite(vals))

    def test_tpms_raw_field(self):
        tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.3)
        raw = tpms.raw_field
        assert callable(raw)
        val = raw(np.array([0.0]), np.array([0.0]), np.array([0.0]))
        assert np.isfinite(val[0])


# ---------------------------------------------------------------------------
# F-rep convenience methods
# ---------------------------------------------------------------------------


class TestTpmsFrepMethods:
    """Test as_sheet, as_upper_skeletal, as_lower_skeletal."""

    def test_as_sheet(self):
        tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.3)
        sheet = tpms.as_sheet(thickness=0.15)
        assert isinstance(sheet, Shape)
        assert sheet.field is not None

    def test_as_sheet_mesh(self):
        tpms = Tpms(
            surface_function=surface_functions.gyroid,
            offset=0.3,
            resolution=20,
        )
        sheet = tpms.as_sheet(thickness=0.15)
        mesh = sheet.generate_surface_mesh(bounds=tpms.bounds, resolution=30)
        assert isinstance(mesh, pv.PolyData)
        assert mesh.n_cells > 0

    def test_as_upper_skeletal(self):
        tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.3)
        skel = tpms.as_upper_skeletal()
        assert isinstance(skel, Shape)
        assert skel.field is not None

    def test_as_lower_skeletal(self):
        tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.3)
        skel = tpms.as_lower_skeletal()
        assert isinstance(skel, Shape)
        assert skel.field is not None


# ---------------------------------------------------------------------------
# Boolean composition
# ---------------------------------------------------------------------------


class TestTpmsBooleans:
    """Test boolean composition of TPMS with other shapes."""

    def test_tpms_and_sphere(self):
        tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.3)
        sphere_func = lambda x, y, z: np.sqrt(x**2 + y**2 + z**2) - 0.4  # noqa: E731
        sphere = from_field(sphere_func, bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5))
        result = tpms & sphere
        assert isinstance(result, Shape)
        mesh = result.generate_surface_mesh(
            bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
            resolution=30,
        )
        assert isinstance(mesh, pv.PolyData)

    def test_sheet_and_sphere(self):
        tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.3)
        sheet = tpms.as_sheet(thickness=0.15)
        sphere_func = lambda x, y, z: np.sqrt(x**2 + y**2 + z**2) - 0.4  # noqa: E731
        sphere = from_field(sphere_func, bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5))
        result = sheet & sphere
        mesh = result.generate_surface_mesh(
            bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
            resolution=30,
        )
        assert isinstance(mesh, pv.PolyData)


# ---------------------------------------------------------------------------
# generate_surface_mesh backward compatibility
# ---------------------------------------------------------------------------


class TestTpmsGenerateVtk:
    """Test that generate_surface_mesh(type_part=...) still works."""

    def test_sheet(self):
        tpms = Tpms(
            surface_function=surface_functions.gyroid,
            offset=0.3,
            resolution=20,
        )
        mesh = tpms.generate_surface_mesh(type_part="sheet")
        assert isinstance(mesh, pv.PolyData)
        assert mesh.n_cells > 0

    def test_surface(self):
        tpms = Tpms(
            surface_function=surface_functions.gyroid,
            offset=0.3,
            resolution=20,
        )
        mesh = tpms.generate_surface_mesh(type_part="surface")
        assert isinstance(mesh, pv.PolyData)
        assert mesh.n_cells > 0

    def test_invalid_type_part(self):
        tpms = Tpms(
            surface_function=surface_functions.gyroid,
            offset=0.3,
        )
        with pytest.raises(ValueError, match="type_part"):
            tpms.generate_surface_mesh(type_part="invalid")

    @pytest.mark.parametrize(
        "surface_fn",
        [surface_functions.gyroid, surface_functions.schwarz_p],
    )
    def test_multiple_surface_functions(self, surface_fn):
        tpms = Tpms(surface_function=surface_fn, offset=0.3, resolution=20)
        mesh = tpms.generate_surface_mesh(type_part="sheet")
        assert mesh.n_cells > 0


# ---------------------------------------------------------------------------
# CylindricalTpms / SphericalTpms
# ---------------------------------------------------------------------------


class TestCurvilinearTpms:
    """Test F-rep field on curvilinear TPMS variants."""

    def test_cylindrical_has_func(self):
        from microgen.shape.tpms import CylindricalTpms

        tpms = CylindricalTpms(
            radius=1.0,
            surface_function=surface_functions.gyroid,
            offset=0.3,
            resolution=10,
        )
        assert tpms.field is not None
        assert tpms.bounds is not None

    def test_cylindrical_evaluate(self):
        from microgen.shape.tpms import CylindricalTpms

        tpms = CylindricalTpms(
            radius=1.0,
            surface_function=surface_functions.gyroid,
            offset=0.3,
            resolution=10,
        )
        x = np.array([1.0, 0.5])
        y = np.array([0.0, 0.5])
        z = np.array([0.0, 0.0])
        vals = tpms.evaluate(x, y, z)
        assert np.all(np.isfinite(vals))

    def test_spherical_has_func(self):
        from microgen.shape.tpms import SphericalTpms

        tpms = SphericalTpms(
            radius=1.0,
            surface_function=surface_functions.gyroid,
            offset=0.3,
            resolution=10,
        )
        assert tpms.field is not None
        assert tpms.bounds is not None

    def test_spherical_evaluate(self):
        from microgen.shape.tpms import SphericalTpms

        tpms = SphericalTpms(
            radius=1.0,
            surface_function=surface_functions.gyroid,
            offset=0.3,
            resolution=10,
        )
        x = np.array([1.0, 0.5])
        y = np.array([0.0, 0.5])
        z = np.array([0.0, 0.5])
        vals = tpms.evaluate(x, y, z)
        assert np.all(np.isfinite(vals))
