"""Tests for the implicit-first paths of :class:`Polyhedron` and lattices.

A convex :class:`Polyhedron` exposes a half-space-intersection SDF
``f(p) = max_k (n_k · (p - c_k))``, suitable for Voronoi cells.
Strut-based lattices expose :meth:`AbstractLattice.generate_implicit`
returning a single composed implicit :class:`Shape`.

Both let downstream :class:`Phase` operate purely on the field —
``Phase.from_shape(...)`` → marching-cubes mesh, grid-quadrature moments,
``pieces`` connectivity — no OCCT required.
"""

import numpy as np

from microgen import Cubic, Phase, Polyhedron
from microgen.shape.shape import Shape

# ruff: noqa: S101


# ------------------------------------------------------------------ Polyhedron


def test_polyhedron_default_tetrahedron_sdf_signs() -> None:
    """SDF is negative at the centroid, ~0 at vertices, positive far outside."""
    poly = Polyhedron()
    f_center = float(
        poly.evaluate(np.array([0.0]), np.array([0.0]), np.array([0.0]))[0]
    )
    f_vertex = float(
        poly.evaluate(np.array([1.0]), np.array([1.0]), np.array([1.0]))[0]
    )
    f_far = float(poly.evaluate(np.array([5.0]), np.array([5.0]), np.array([5.0]))[0])
    assert f_center < 0
    assert abs(f_vertex) < 1e-6
    assert f_far > 0


def test_polyhedron_volume_via_phase_matches_analytic() -> None:
    """Tetrahedron with vertices (±1,±1,±1)-pattern has volume 8/3."""
    poly = Polyhedron()
    phase = Phase.from_shape(poly, resolution=80)
    total_vol = sum(p.volume for p in phase.pieces)
    assert np.isclose(total_vol, 8.0 / 3.0, rtol=0.05)


def test_polyhedron_phase_is_one_piece() -> None:
    """A single convex polyhedron yields one connected component."""
    poly = Polyhedron()
    phase = Phase.from_shape(poly, resolution=50)
    assert len(phase.pieces) == 1


def test_polyhedron_translates_correctly() -> None:
    """When ``center`` is non-zero, the field's COM matches the center."""
    poly = Polyhedron(center=(2.0, -1.0, 0.5))
    phase = Phase.from_shape(poly, resolution=50)
    assert np.allclose(phase.center_of_mass, (2.0, -1.0, 0.5), atol=0.05)


# ------------------------------------------------------------------ Lattice


def test_cubic_lattice_implicit_returns_shape() -> None:
    """``generate_implicit()`` returns a composed :class:`Shape`."""
    lat = Cubic(strut_radius=0.05)
    implicit = lat.generate_implicit()
    assert isinstance(implicit, Shape)
    assert implicit.func is not None
    assert implicit.bounds is not None


def test_cubic_lattice_implicit_one_connected_piece() -> None:
    """The 12 edge struts share corner endpoints — one connected piece."""
    lat = Cubic(strut_radius=0.05)
    phase = Phase.from_shape(lat.generate_implicit(), resolution=60)
    assert len(phase.pieces) == 1


def test_cubic_lattice_implicit_volume_approaches_cad() -> None:
    """Grid quadrature on the implicit lattice converges toward the CAD volume."""
    lat = Cubic(strut_radius=0.05)
    phase = Phase.from_shape(lat.generate_implicit(), resolution=100)
    grid_vol = sum(p.volume for p in phase.pieces)
    cad_vol = lat.cad_shape.volume()
    # 100³ quadrature on a 0.05-radius lattice: expect ~5% error.
    assert np.isclose(grid_vol, cad_vol, rtol=0.10)


def test_cubic_lattice_implicit_center_of_mass_at_origin() -> None:
    """Symmetric cubic cell — COM should be near the cell center."""
    lat = Cubic(strut_radius=0.05)
    phase = Phase.from_shape(lat.generate_implicit(), resolution=60)
    assert np.allclose(phase.center_of_mass, (0.0, 0.0, 0.0), atol=1e-2)
