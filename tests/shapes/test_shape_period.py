"""Tests for ``Shape.period``: the intrinsic-periodicity attribute.

When set, ``shape.period == (Lx, Ly, Lz)`` is a data-structure invariant
guaranteeing ``shape.evaluate(p + L) == shape.evaluate(p)`` along each axis.
``Tpms`` / ``Spinodoid`` populate it from ``cell_size * repeat_cell`` /
``cell_size`` respectively; free shapes built via ``from_field`` or boolean
composition leave it ``None``.
"""

# ruff: noqa: S101

from __future__ import annotations

import numpy as np

from microgen import Box, Sphere, Spinodoid, Tpms, surface_functions
from microgen.shape.implicit_ops import from_field
from microgen.shape.shape import Shape


def test_bare_shape_period_is_none() -> None:
    """A free ``Shape(func=...)`` has no intrinsic period."""
    s = Shape(func=lambda x, y, z: x * x + y * y + z * z - 1.0)
    assert s.period is None


def test_box_sphere_have_no_period() -> None:
    """Non-periodic primitives expose ``period is None``."""
    assert Box().period is None
    assert Sphere().period is None


def test_tpms_period_matches_cell_size() -> None:
    """``Tpms.period`` equals ``cell_size`` (per-axis)."""
    tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.5, cell_size=2.0)
    assert tpms.period == (2.0, 2.0, 2.0)


def test_tpms_anisotropic_cell_size_period() -> None:
    """``Tpms.period`` is per-axis when ``cell_size`` is a tuple."""
    tpms = Tpms(
        surface_function=surface_functions.gyroid,
        offset=0.5,
        cell_size=(1.0, 2.0, 3.0),
    )
    assert tpms.period == (1.0, 2.0, 3.0)


def test_spinodoid_period_matches_cell_size() -> None:
    """``Spinodoid.period`` equals ``cell_size``."""
    sp = Spinodoid(offset=0.0, cell_size=1.5, resolution=8, seed=42)
    assert sp.period == (1.5, 1.5, 1.5)


def test_from_field_propagates_no_period() -> None:
    """``from_field`` produces a free shape with ``period is None``."""
    s = from_field(lambda x, y, z: x + y + z)
    assert s.period is None


def test_boolean_composition_does_not_inherit_period() -> None:
    """``a | b`` produces a free shape; intrinsic periodicity is not preserved
    by composition (the union of two periodic fields may or may not be
    periodic — we don't claim it is).
    """
    tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.5)
    sphere = Sphere(radius=0.5)
    combined = tpms | sphere
    assert combined.period is None


def test_tpms_field_actually_periodic_at_declared_period() -> None:
    """If ``shape.period == (Lx, Ly, Lz)`` then ``f(p + L) == f(p)``."""
    tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.5, cell_size=1.0)
    lx, ly, lz = tpms.period
    rng = np.random.default_rng(seed=0)
    pts = rng.uniform(-0.5, 0.5, size=(20, 3))
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    f0 = tpms.evaluate(x, y, z)
    fp = tpms.evaluate(x + lx, y + ly, z + lz)
    # TPMS field is approximately periodic after SDF normalization (not bit-exact);
    # tight but not exact tolerance.
    assert np.allclose(f0, fp, atol=1e-6)


def test_spinodoid_field_periodic_at_declared_period() -> None:
    """Spinodoid's reciprocal-lattice modes give (float-precision) periodicity.

    Bit-exactness is asserted by ``test_spinodoid_field_is_bit_exact_periodic``
    in tests/test_spinodoid.py on the raw ``_frep.evaluate``; via ``Shape.evaluate``
    the sign-flip + shifted-coord path picks up ULP-level drift (~1e-15).
    """
    sp = Spinodoid(offset=0.0, cell_size=1.0, resolution=8, seed=42)
    lx, ly, lz = sp.period
    rng = np.random.default_rng(seed=1)
    pts = rng.uniform(0.0, 1.0, size=(20, 3))
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    f0 = sp.evaluate(x, y, z)
    fp = sp.evaluate(x + lx, y + ly, z + lz)
    assert np.allclose(f0, fp, atol=1e-12, rtol=0)
