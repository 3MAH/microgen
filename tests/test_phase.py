"""Test for the phase module."""

import numpy as np

from microgen import Phase, Rve, rasterPhase
from microgen.shape import Box, Ellipsoid, Sphere

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/


def test_phase_sphere_rasterize_must_have_correct_number_of_solids() -> None:
    """Rasterize Sphere by 3x3x3 grid must have 27 solids."""
    rve = Rve(dim=1)
    sphere = Sphere(radius=0.5).generate()
    grid = [3 for _ in range(3)]

    phase = Phase(shape=sphere)
    phases = phase.rasterize(rve=rve, grid=grid)
    raster = rasterPhase(phase=phase, rve=rve, grid=grid, phasePerRaster=False)

    assert len(phases) == np.prod(grid)
    assert len(phases) == len(raster.solids)


def test_phase_rasterize_phase_per_raster_should_return_the_right_object() -> None:
    """Test the Phase class with a shape."""
    rve = Rve(dim=1)
    sphere = Sphere(radius=0.5).generate()
    phase = Phase(shape=sphere)

    grid = [3, 3, 3]
    phases = phase.rasterize(rve, grid, phasePerRaster=True)
    assert isinstance(phases, list)
    assert len(phases) == np.prod(grid)

    grid = [2, 2, 2]
    phase.rasterize(rve, grid, phasePerRaster=False)
    assert len(phase.solids) == np.prod(grid)


def test_phase_repeat_should_repeat_the_shape_in_the_rve() -> None:
    """Test the Phase class with a shape."""
    rve = Rve(dim=1)
    box = Box(dim=(1.0, 1.0, 1.0)).generate()
    volume_before = box.Volume()
    phase = Phase(shape=box)

    repeat = (1, 2, 1)
    phase.repeat(rve, grid=repeat)
    assert len(phase.solids) == np.prod(repeat)
    assert phase.shape.Volume() == 2.0 * volume_before


def test_phase_rescale_should_change_the_size_of_the_shape() -> None:
    """Test the Phase class with a shape."""
    radius = 1.0
    scale = 1.5
    phase = Phase(shape=Sphere(radius=radius).generate())
    volume_before = phase.shape.Volume()
    phase.rescale(scale)
    assert np.isclose(phase.shape.Volume(), volume_before * scale**3, rtol=1e-2)


def test_phase_translate_should_shift_centers_corresponding_to_the_translation() -> (
    None
):
    """Test the Phase class with a shape."""
    center = (1.0, 0.5, -0.5)
    phase = Phase(shape=Ellipsoid(center=center, radii=(0.15, 0.31, 0.4)).generate())
    phase.translate((1, 0, 0))
    assert np.allclose(phase.centerOfMass, (2.0, 0.5, -0.5), rtol=1e-4)
    phase.translate(np.array([0, 1, 1]))
    assert np.allclose(phase.centerOfMass, (2.0, 1.5, 0.5), rtol=1e-4)


def test_phase_center_of_mass_should_return_the_right_values() -> None:
    """Test the Phase class with a shape."""
    center = (1.0, 0.5, -0.5)
    phase = Phase(shape=Ellipsoid(center=center, radii=(0.15, 0.31, 0.4)).generate())
    assert np.allclose(phase.centerOfMass, center, rtol=1e-4)
    assert np.allclose(
        phase.getCenterOfMass(compute=False),
        phase.centerOfMass,
        rtol=1e-4,
    )


def test_phase_inertia_matrix_should_return_the_right_values() -> None:
    """Test the Phase class with a shape."""
    radius = 1.5
    phase = Phase(shape=Sphere(radius=radius).generate())
    assert np.allclose(
        phase.inertiaMatrix,
        np.diag([phase.shape.Volume() * 2 / 5 * radius**2] * 3),
    )
    assert np.allclose(phase.getInertiaMatrix(compute=False), phase.inertiaMatrix)


def test_phase_solids_and_shape_properties_should_return_the_right_values() -> None:
    """Test the Phase class with a shape."""
    ellipsoid = Ellipsoid(radii=(0.15, 0.31, 0.4)).generate()
    phase = Phase(shape=ellipsoid)
    assert phase.solids == [ellipsoid]
    assert phase.shape == ellipsoid


def test_phase_empty_should_have_empty_shape_and_solids() -> None:
    """Test the Phase class with no shape."""
    void_phase = Phase()
    assert void_phase.shape is None
    assert void_phase.solids == []
