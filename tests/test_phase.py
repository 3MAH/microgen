"""Tests for the Phase 2.0 module.

Phase is now field-first (implicit-first) with CAD as an optional view.
Construction goes through ``Phase.from_cad`` for CAD-backed phases and
``Phase.from_shape`` for field-backed ones.  Transforms are immutable —
``translated`` / ``scaled`` / ``tiled`` return new :class:`Phase`
instances rather than mutating in place.
"""

import numpy as np

from microgen import Phase, Rve, raster_phase
from microgen.shape import Box, Ellipsoid, Sphere

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/


def test_phase_sphere_rasterize_must_have_correct_number_of_solids() -> None:
    """Rasterize Sphere by 3x3x3 grid must have 27 phases / 27 sub-solids."""
    rve = Rve(dim=1)
    sphere = Sphere(radius=0.5).generate_cad()
    grid = [3 for _ in range(3)]

    phase = Phase.from_cad(sphere)
    phases = raster_phase(phase=phase, rve=rve, grid=grid)
    raster = raster_phase(phase=phase, rve=rve, grid=grid, phase_per_raster=False)

    assert len(phases) == np.prod(grid)
    assert len(raster.cad.solids()) == np.prod(grid)


def test_phase_rasterize_phase_per_raster_should_return_the_right_object() -> None:
    """Test rastering via the operations free function."""
    rve = Rve(dim=1)
    sphere = Sphere(radius=0.5).generate_cad()
    phase = Phase.from_cad(sphere)

    grid = [3, 3, 3]
    phases = raster_phase(phase=phase, rve=rve, grid=grid, phase_per_raster=True)
    assert isinstance(phases, list)
    assert len(phases) == np.prod(grid)

    grid = [2, 2, 2]
    rastered = raster_phase(phase=phase, rve=rve, grid=grid, phase_per_raster=False)
    assert len(rastered.cad.solids()) == np.prod(grid)


def test_phase_tiled_should_repeat_the_shape_in_the_rve() -> None:
    """Phase.tiled returns a new Phase with the geometry tiled on a grid."""
    rve = Rve(dim=1)
    box = Box(dim=(1.0, 1.0, 1.0)).generate_cad()
    volume_before = box.volume()
    phase = Phase.from_cad(box)

    repeat = (1, 2, 1)
    tiled = phase.tiled(rve, grid=repeat)
    assert len(tiled.cad.solids()) == np.prod(repeat)
    assert np.isclose(tiled.cad.volume(), 2.0 * volume_before)


def test_phase_scaled_should_change_the_size_of_the_shape() -> None:
    """Phase.scaled returns a new Phase with cubed volume."""
    radius = 1.0
    scale = 1.5
    phase = Phase.from_cad(Sphere(radius=radius).generate_cad())
    volume_before = phase.cad.volume()
    scaled = phase.scaled(scale)
    assert np.isclose(scaled.cad.volume(), volume_before * scale**3, rtol=1e-2)


def test_phase_translated_should_shift_centers() -> None:
    """Phase.translated returns a new Phase with the COM shifted."""
    center = (1.0, 0.5, -0.5)
    phase = Phase.from_cad(
        Ellipsoid(center=center, radii=(0.15, 0.31, 0.4)).generate_cad()
    )
    moved = phase.translated((1, 0, 0))
    assert np.allclose(moved.center_of_mass, (2.0, 0.5, -0.5), rtol=1e-4)
    moved2 = moved.translated(np.array([0, 1, 1]))
    assert np.allclose(moved2.center_of_mass, (2.0, 1.5, 0.5), rtol=1e-4)


def test_phase_center_of_mass_should_return_the_right_values() -> None:
    """``phase.center_of_mass`` is the BREP volumetric COM."""
    center = (1.0, 0.5, -0.5)
    phase = Phase.from_cad(
        Ellipsoid(center=center, radii=(0.15, 0.31, 0.4)).generate_cad()
    )
    assert np.allclose(phase.center_of_mass, center, rtol=1e-4)
    # Cached: accessing twice yields the same array object.
    assert phase.center_of_mass is phase.center_of_mass


def test_phase_inertia_matrix_should_return_the_right_values() -> None:
    """Inertia tensor of a uniform sphere is ``(2/5) m R²`` on the diagonal."""
    radius = 1.5
    phase = Phase.from_cad(Sphere(radius=radius).generate_cad())
    assert np.allclose(
        phase.inertia_matrix,
        np.diag([phase.cad.volume() * 2 / 5 * radius**2] * 3),
    )
    # Cached.
    assert phase.inertia_matrix is phase.inertia_matrix


def test_phase_cad_and_pieces_properties() -> None:
    """``phase.cad`` exposes the CAD shape; ``phase.pieces`` enumerates solids."""
    ellipsoid = Ellipsoid(radii=(0.15, 0.31, 0.4)).generate_cad()
    phase = Phase.from_cad(ellipsoid)
    assert phase.cad is ellipsoid
    pieces = phase.pieces
    assert len(pieces) == 1
    # The Piece carries its CAD payload.
    assert pieces[0].cad is not None
    assert pieces[0].volume > 0.0


def test_phase_empty_is_reported_as_such() -> None:
    """``Phase()`` (no field/cad/mesh) is empty."""
    void_phase = Phase()
    assert void_phase.is_empty
    assert void_phase.field is None
    assert void_phase.bounds is None


def test_phase_from_shape_field_backed() -> None:
    """``Phase.from_shape`` builds a field-backed phase from an implicit Shape."""
    sphere = Sphere(radius=0.5)
    phase = Phase.from_shape(sphere)
    assert phase.field is not None
    assert phase.bounds is not None
    # Center of mass via grid quadrature on a centered sphere ≈ origin.
    assert np.allclose(phase.center_of_mass, (0.0, 0.0, 0.0), atol=2e-2)
    # The sphere is one connected piece.
    assert len(phase.pieces) == 1


def test_phase_from_shape_propagates_period() -> None:
    """When the source Shape has a period, the Phase carries it."""
    from microgen import Tpms, surface_functions

    tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.5, cell_size=1.0)
    phase = Phase.from_shape(tpms)
    assert phase.period == (1.0, 1.0, 1.0)
