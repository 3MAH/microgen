"""Test for the phase module."""

import numpy as np

from microgen import Phase, Rve, raster_phase
from microgen.shape import Box, Ellipsoid, Sphere

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/


def test_phase_sphere_rasterize_must_have_correct_number_of_solids() -> None:
    """Rasterize Sphere by 3x3x3 grid must have 27 solids."""
    rve = Rve(dim=1)
    sphere = Sphere(radius=0.5).generate_cad()
    grid = [3 for _ in range(3)]

    phase = Phase(shape=sphere)
    phases = phase.rasterize(rve=rve, grid=grid)
    raster = raster_phase(phase=phase, rve=rve, grid=grid, phase_per_raster=False)

    assert len(phases) == np.prod(grid)
    assert len(phases) == len(raster.solids)


def test_phase_rasterize_phase_per_raster_should_return_the_right_object() -> None:
    """Test the Phase class with a shape."""
    rve = Rve(dim=1)
    sphere = Sphere(radius=0.5).generate_cad()
    phase = Phase(shape=sphere)

    grid = [3, 3, 3]
    phases = phase.rasterize(rve, grid, phase_per_raster=True)
    assert isinstance(phases, list)
    assert len(phases) == np.prod(grid)

    grid = [2, 2, 2]
    phase.rasterize(rve, grid, phase_per_raster=False)
    assert len(phase.solids) == np.prod(grid)


def test_phase_repeat_should_repeat_the_shape_in_the_rve() -> None:
    """Test the Phase class with a shape."""
    rve = Rve(dim=1)
    box = Box(dim=(1.0, 1.0, 1.0)).generate_cad()
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
    phase = Phase(shape=Sphere(radius=radius).generate_cad())
    volume_before = phase.shape.Volume()
    phase.rescale(scale)
    assert np.isclose(phase.shape.Volume(), volume_before * scale**3, rtol=1e-2)


def test_phase_translate_should_shift_centers_corresponding_to_the_translation() -> (
    None
):
    """Test the Phase class with a shape."""
    center = (1.0, 0.5, -0.5)
    phase = Phase(
        shape=Ellipsoid(center=center, radii=(0.15, 0.31, 0.4)).generate_cad()
    )
    phase.translate((1, 0, 0))
    assert np.allclose(phase.center_of_mass, (2.0, 0.5, -0.5), rtol=1e-4)
    phase.translate(np.array([0, 1, 1]))
    assert np.allclose(phase.center_of_mass, (2.0, 1.5, 0.5), rtol=1e-4)


def test_phase_center_of_mass_should_return_the_right_values() -> None:
    """Test the Phase class with a shape."""
    center = (1.0, 0.5, -0.5)
    phase = Phase(
        shape=Ellipsoid(center=center, radii=(0.15, 0.31, 0.4)).generate_cad()
    )
    assert np.allclose(phase.center_of_mass, center, rtol=1e-4)
    assert np.allclose(
        phase.get_center_of_mass(compute=False),
        phase.center_of_mass,
        rtol=1e-4,
    )


def test_phase_inertia_matrix_should_return_the_right_values() -> None:
    """Test the Phase class with a shape."""
    radius = 1.5
    phase = Phase(shape=Sphere(radius=radius).generate_cad())
    assert np.allclose(
        phase.inertia_matrix,
        np.diag([phase.shape.Volume() * 2 / 5 * radius**2] * 3),
    )
    assert np.allclose(phase.get_inertia_matrix(compute=False), phase.inertia_matrix)


def test_phase_solids_and_shape_properties_should_return_the_right_values() -> None:
    """Test the Phase class with a shape."""
    ellipsoid = Ellipsoid(radii=(0.15, 0.31, 0.4)).generate_cad()
    phase = Phase(shape=ellipsoid)
    # Phase.solids returns raw TopoDS_Solid objects; identity-by-`==` doesn't
    # apply, so check topological sameness against the wrapped solid instead.
    assert len(phase.solids) == 1
    assert phase.solids[0].IsSame(ellipsoid.wrapped)
    assert phase.shape is ellipsoid


def test_phase_empty_should_have_empty_shape_and_solids() -> None:
    """Test the Phase class with no shape."""
    void_phase = Phase()
    assert void_phase.shape is None
    assert void_phase.solids == []
