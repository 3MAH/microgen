"""Tests for the base-class ``Shape.generate_surface_mesh`` / ``generate_volume_mesh``.

Primitives (``Sphere``, ``Box``, ``Capsule``, ``Cylinder``, ``Ellipsoid``)
override ``generate_surface_mesh`` with native pyvista renderers, but
**inherit** ``generate_volume_mesh`` from the base ``Shape`` (the F-rep
sample + clip path). This file covers:

- the inherited primitive path lands at the declared ``center`` exactly
  (regression net for the prior double-translate bug);
- the same path works for anonymous shapes built via
  :func:`microgen.shape.implicit_ops.from_field` and via boolean
  composition;
- the grid sample is cached across ``generate_surface_mesh`` +
  ``generate_volume_mesh`` calls on the same instance;
- ``center`` / ``orientation`` are read-only on the base, but propagate
  correctly through ``translate`` / ``rotate`` / ``scale``.
"""

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/

from __future__ import annotations

import numpy as np
import pytest
import pyvista as pv
from scipy.spatial.transform import Rotation

from microgen import Box, Capsule, Cylinder, Ellipsoid, Sphere
from microgen.shape.implicit_ops import from_field
from microgen.shape.shape import Shape


# ---------------------------------------------------------------------------
# Default volume-mesh path on primitives (regression for center double-apply)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("ctor", "kwargs"),
    [
        (Sphere, {"radius": 1.0, "center": (5, 0, 0)}),
        (Box, {"dim": (1, 1, 1), "center": (0, 3, 0)}),
        (Capsule, {"radius": 0.3, "height": 1.0, "center": (2, 2, 0)}),
        (Cylinder, {"radius": 0.5, "height": 1.0, "center": (-3, 0, 1)}),
        (Ellipsoid, {"radii": (1.0, 0.5, 0.25), "center": (0, 0, 4)}),
    ],
)
def test_primitive_volume_mesh_lands_at_center(ctor, kwargs):
    """Inherited ``generate_volume_mesh`` lands at the declared ``center``."""
    s = ctor(**kwargs)
    mesh = s.generate_volume_mesh(resolution=20)
    assert mesh.n_points > 0
    actual_center = np.asarray(mesh.center, dtype=np.float64)
    expected_center = np.asarray(kwargs["center"], dtype=np.float64)
    # Voxel-grid centers within half a grid cell.
    assert np.linalg.norm(actual_center - expected_center) < 0.5


def test_primitive_volume_mesh_is_unstructured_grid():
    """The default path returns ``pv.UnstructuredGrid``, not ``PolyData``."""
    mesh = Sphere(radius=1.0).generate_volume_mesh(resolution=20)
    assert isinstance(mesh, pv.UnstructuredGrid)
    assert mesh.n_cells > 0


# ---------------------------------------------------------------------------
# Anonymous shapes via from_field and boolean composition
# ---------------------------------------------------------------------------


def test_from_field_volume_mesh():
    """A shape built via ``from_field`` exposes the default volume-mesh path."""
    s = from_field(
        field=lambda x, y, z: np.sqrt(x**2 + y**2 + z**2) - 1.0,
        bounds=(-1.2, 1.2, -1.2, 1.2, -1.2, 1.2),
    )
    mesh = s.generate_volume_mesh(resolution=20)
    assert isinstance(mesh, pv.UnstructuredGrid)
    assert mesh.n_cells > 0


def test_composed_shape_volume_mesh():
    """Boolean composition produces a Shape with a working default path."""
    a = Sphere(radius=0.6, center=(-0.4, 0, 0))
    b = Sphere(radius=0.6, center=(0.4, 0, 0))
    union = a | b
    mesh = union.generate_volume_mesh(resolution=30)
    assert mesh.n_points > 0


# ---------------------------------------------------------------------------
# Grid-cache: consecutive surface + volume calls reuse one sampling
# ---------------------------------------------------------------------------


def test_grid_cache_avoids_resampling():
    """Two calls on the same instance hit the cache the second time."""
    s = from_field(
        field=lambda x, y, z: x**2 + y**2 + z**2 - 1.0,
        bounds=(-1.2, 1.2, -1.2, 1.2, -1.2, 1.2),
    )
    s.generate_surface_mesh(resolution=20)
    # Internal cache populated with the (bounds, resolution) key.
    assert s._grid_cache  # noqa: SLF001
    cached = next(iter(s._grid_cache.values()))  # noqa: SLF001
    # Volume call must reuse the same grid object.
    s.generate_volume_mesh(resolution=20)
    assert next(iter(s._grid_cache.values())) is cached  # noqa: SLF001


def test_grid_cache_keyed_on_resolution():
    """A different ``resolution`` produces a fresh cache entry."""
    s = from_field(
        field=lambda x, y, z: x**2 + y**2 + z**2 - 1.0,
        bounds=(-1.2, 1.2, -1.2, 1.2, -1.2, 1.2),
    )
    s.generate_surface_mesh(resolution=20)
    s.generate_surface_mesh(resolution=40)
    assert len(s._grid_cache) == 2  # noqa: SLF001


# ---------------------------------------------------------------------------
# Read-only center / orientation on the base
# ---------------------------------------------------------------------------


def test_center_is_read_only_on_shape():
    """Setting ``.center`` on a bare ``Shape`` raises ``AttributeError``."""
    s = from_field(lambda x, y, z: x, bounds=(0, 1, 0, 1, 0, 1))
    with pytest.raises(AttributeError):
        s.center = (1, 2, 3)


def test_orientation_is_read_only_on_shape():
    """Setting ``.orientation`` on a bare ``Shape`` raises ``AttributeError``."""
    s = from_field(lambda x, y, z: x, bounds=(0, 1, 0, 1, 0, 1))
    with pytest.raises(AttributeError):
        s.orientation = Rotation.identity()


# ---------------------------------------------------------------------------
# translate / rotate / scale propagate center + orientation
# ---------------------------------------------------------------------------


def test_translated_propagates_center():
    """``translated`` shifts the reported center by the offset."""
    s = Sphere(radius=1.0, center=(1.0, 2.0, 3.0))
    shifted = s.translated((4.0, -1.0, 2.0))
    assert tuple(shifted.center) == pytest.approx((5.0, 1.0, 5.0))


def test_translated_preserves_orientation():
    """``translated`` does not touch orientation."""
    s = Sphere(radius=1.0, orientation=Rotation.from_euler("z", 30, degrees=True))
    shifted = s.translated((1.0, 0.0, 0.0))
    np.testing.assert_allclose(
        shifted.orientation.as_matrix(),
        s.orientation.as_matrix(),
    )


def test_rotated_propagates_center_and_orientation():
    """``rotated`` rotates the reported center about the world origin."""
    s = Sphere(radius=1.0, center=(1.0, 0.0, 0.0))
    rotated = s.rotated((0.0, 0.0, 90.0), convention="xyz")
    # Rotation of (1, 0, 0) by 90° about z lands at (0, 1, 0).
    np.testing.assert_allclose(rotated.center, (0.0, 1.0, 0.0), atol=1e-9)


def test_scaled_propagates_center():
    """``scaled`` scales the reported center by the same factor."""
    s = Sphere(radius=1.0, center=(2.0, 0.0, 0.0))
    scaled = s.scaled(3.0)
    assert tuple(scaled.center) == pytest.approx((6.0, 0.0, 0.0))


# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------


def test_volume_mesh_without_func_raises():
    """Bare ``Shape()`` (no ``_func``) raises a useful error from volume_mesh."""
    s = Shape()
    with pytest.raises(NotImplementedError, match="No implicit field"):
        s.generate_volume_mesh(bounds=(0, 1, 0, 1, 0, 1))


def test_volume_mesh_without_bounds_raises():
    """A shape with a func but no bounds raises when bounds aren't passed."""
    s = Shape(field=lambda x, y, z: x**2 + y**2 + z**2 - 1.0)
    with pytest.raises(ValueError, match="Bounds must be provided"):
        s.generate_volume_mesh()
