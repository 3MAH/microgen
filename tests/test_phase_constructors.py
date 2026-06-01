"""Tests for the new Tier 2 :class:`Phase` constructors.

- :meth:`Phase.from_implicit` — sugar over the field-first constructor that
  derives ``bounds`` from a :class:`Rve`.
- :meth:`Phase.from_grid` — wraps a pre-sampled
  :class:`pyvista.StructuredGrid` (useful for expensive-to-evaluate fields).
- :meth:`Phase.rotated` — the missing immutable transform paired with
  ``translated`` / ``scaled`` / ``tiled``.
"""

import numpy as np
import pyvista as pv

from microgen import Box, Phase, Rve, Sphere

# ruff: noqa: S101


def _sphere_sdf(x, y, z, r=0.5):
    return np.sqrt(x * x + y * y + z * z) - r


def test_from_implicit_derives_bounds_from_rve() -> None:
    """``Phase.from_implicit(f, rve)`` derives bounds from ``rve.min/max_point``."""
    rve = Rve(center=(0.0, 0.0, 0.0), dim=2.0)
    phase = Phase.from_implicit(_sphere_sdf, rve, resolution=30)
    assert phase.bounds == (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
    assert phase.field is not None
    assert phase.period is None
    # Centred sphere → COM at origin.
    assert np.allclose(phase.center_of_mass, (0.0, 0.0, 0.0), atol=0.05)


def test_from_implicit_propagates_period_and_name() -> None:
    """Period and name kwargs flow through unchanged."""
    rve = Rve(dim=1.0)
    phase = Phase.from_implicit(
        _sphere_sdf, rve, period=(1.0, 1.0, 1.0), name="my_phase"
    )
    assert phase.period == (1.0, 1.0, 1.0)
    assert phase.name == "my_phase"


def test_from_grid_returns_input_grid_untouched() -> None:
    """``Phase.from_grid`` reuses the provided grid as its source of truth."""
    nx = 30
    xi = np.linspace(-1, 1, nx)
    x, y, z = np.meshgrid(xi, xi, xi, indexing="ij")
    sg = pv.StructuredGrid(x, y, z)
    sg["implicit"] = _sphere_sdf(x, y, z).ravel(order="F")

    phase = Phase.from_grid(sg, scalars="implicit")
    # Cached grid is the same object (no resampling).
    assert phase.grid() is sg
    assert phase.bounds == (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
    # COM of the centred sphere ≈ origin.
    assert np.allclose(phase.center_of_mass, (0.0, 0.0, 0.0), atol=0.05)


def test_from_grid_rejects_missing_scalar() -> None:
    """Asking for a scalar name that isn't on the grid raises ``ValueError``."""
    import pytest

    sg = pv.StructuredGrid(
        *np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij")
    )
    sg["other"] = np.zeros(sg.n_points)
    with pytest.raises(ValueError, match="no point scalar named 'implicit'"):
        Phase.from_grid(sg)


def test_phase_rotated_field_backed_swaps_axes() -> None:
    """A box along x rotated 90° about z lands along y (field-backed)."""
    box = Box(dim=(2.0, 0.4, 0.4))
    phase = Phase.from_shape(box, resolution=40)
    # f(0.5, 0, 0) < 0 before, > 0 after; f(0, 0.5, 0) < 0 after.
    f_before = phase.field(np.array([0.5]), np.array([0.0]), np.array([0.0]))[0]
    rotated = phase.rotated((0, 0, 90), convention="xyz")
    f_after = rotated.field(np.array([0.5]), np.array([0.0]), np.array([0.0]))[0]
    f_y = rotated.field(np.array([0.0]), np.array([0.5]), np.array([0.0]))[0]
    assert f_before < 0
    assert f_after > 0
    assert f_y < 0


def test_phase_rotated_field_backed_invalidates_period() -> None:
    """Rotation breaks axis-aligned periodicity; ``period`` resets to ``None``."""
    from microgen import Tpms, surface_functions

    tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.3, cell_size=1.0)
    phase = Phase.from_shape(tpms)
    assert phase.period == (1.0, 1.0, 1.0)
    rotated = phase.rotated((30, 0, 0))
    assert rotated.period is None


def test_phase_rotated_cad_backed_preserves_volume() -> None:
    """CAD-backed rotation produces a Phase with the same volume."""
    sph = Phase.from_cad(Sphere(radius=0.5).generate_cad())
    rotated = sph.rotated((0, 0, 45))
    # CAD-backed: same volume (rigid transform).
    assert np.isclose(rotated.cad.volume(), sph.cad.volume(), rtol=1e-6)


def test_phase_rotated_empty_raises() -> None:
    """An empty Phase cannot be rotated."""
    import pytest

    with pytest.raises(ValueError, match="Cannot rotate an empty Phase"):
        Phase().rotated((10, 0, 0))
