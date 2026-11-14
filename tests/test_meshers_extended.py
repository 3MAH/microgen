"""Field and domain controls for sampled, mapped and clipped TPMS."""

import numpy as np
import pyvista as pv
import pytest

from microgen import CylindricalTpms, Infill, NormedDistance, SphericalTpms, Sweep, Tpms
from microgen.shape.surface_functions import gyroid


@pytest.mark.parametrize("kind", [CylindricalTpms, SphericalTpms])
def test_mapped_raw_field_matches_original_parametric_grid(kind):
    shape = kind(radius=2, surface_function=gyroid, offset=1.1, resolution=8)
    np.testing.assert_allclose(
        shape.raw_field(*shape.grid.points.T), shape.grid["surface"], atol=1e-12
    )


def test_sampled_offset_is_frozen_and_batch_independent():
    shape = Tpms(gyroid, offset=1.1, resolution=8)
    samples = 1.1 + shape.grid.points @ np.array([0.1, 0.2, 0.3])
    shape.offset = samples
    points = np.random.default_rng(4).uniform(-0.5, 0.5, (25, 3))
    expected = 1.1 + points @ np.array([0.1, 0.2, 0.3])
    assert shape._offset_func is not None
    np.testing.assert_allclose(shape._offset_func(*points.T), expected, atol=1e-12)
    samples[:] = 99
    np.testing.assert_allclose(shape._offset_func(*points[:1].T), expected[:1])
    assert shape.offset.max() < 2


def test_grading_uses_original_reference_grid_normalization():
    grading = NormedDistance(
        pv.Plane(direction=(1, 0, 0), i_size=4, j_size=4),
        boundary_offset=1.2,
        furthest_offset=0.8,
    )
    shape = Tpms(gyroid, offset=grading, resolution=8)
    expected = grading.compute_offset(shape.grid)
    np.testing.assert_array_equal(shape.offset, expected)
    points = shape.grid.points[[0, 9, 20]]
    assert shape._offset_func is not None
    batch = shape._offset_func(*points.T)
    alone = shape._offset_func(*points[:1].T)
    np.testing.assert_allclose(batch, expected[[0, 9, 20]], atol=1e-12)
    np.testing.assert_array_equal(alone, batch[:1])


@pytest.mark.parametrize("samples", [np.ones(3), np.full(8**3, np.nan)])
def test_invalid_sampled_offset_is_rejected(samples):
    shape = Tpms(gyroid, offset=1.1, resolution=8)
    with pytest.raises(ValueError, match="sample"):
        shape.offset = samples


@pytest.mark.parametrize("kind", [CylindricalTpms, SphericalTpms])
def test_mapped_mesh_has_positive_tetrahedra_and_retains_parameters(kind):
    pytest.importorskip("meshers")
    shape = kind(radius=2, surface_function=gyroid, offset=1.1, resolution=8)
    mesh = shape.generate_volume_mesh(backend="meshers", cells=16)
    assert mesh.n_cells > 0
    corners = mesh.points[mesh.cells_dict[pv.CellType.TETRA]]
    assert np.linalg.det(corners[:, 1:] - corners[:, :1]).min() > 0
    assert "meshers_parameters" in mesh.point_data
    assert mesh.extract_surface(algorithm=None).n_open_edges == 0


def test_full_cylindrical_seam_is_welded():
    pytest.importorskip("meshers")
    shape = CylindricalTpms(
        radius=2,
        surface_function=gyroid,
        density=1,
        repeat_cell=(1, 0, 1),
        resolution=8,
    )
    mesh = shape.generate_volume_mesh(backend="meshers", cells=(8, 64, 8))
    assert mesh.extract_surface(algorithm=None).n_open_edges == 0
    assert "periodic_pairs_y" not in mesh.field_data
    assert mesh.volume == pytest.approx(4 * np.pi, rel=0.01)


def test_infill_retains_its_actual_envelope():
    pytest.importorskip("meshers")
    obj = pv.Sphere(radius=0.45, theta_resolution=32, phi_resolution=32)
    shape = Infill(
        obj=obj, surface_function=gyroid, offset=1.1, repeat_cell=1, resolution=12
    )
    mesh = shape.generate_volume_mesh(backend="meshers", cells=24)
    assert 0 < mesh.volume < obj.volume
    assert np.linalg.norm(mesh.points, axis=1).max() <= 0.451
    assert "envelope" in mesh.field_data["meshers_constraint_names"]


def test_annular_sweep_map_and_field_agree():
    shape = Sweep(
        curve_points=np.array([[0.0, 0, 0], [0.0, 0, 0.5], [0.0, 0, 1]]),
        surface_function=gyroid,
        radial_max=0.4,
        radial_min=0.2,
        offset=1.1,
        resolution=8,
    )
    np.testing.assert_allclose(
        shape.raw_field(*shape.grid.points.T), shape.grid["surface"], atol=1e-10
    )
    pytest.importorskip("meshers")
    # Isolate the regular sweep domain from TPMS cut-cell quality limitations.
    domain = Sweep(
        curve_points=shape.curve,
        surface_function=gyroid,
        radial_max=0.4,
        radial_min=0.2,
        density=1,
        resolution=8,
    )
    mesh = domain.generate_volume_mesh(backend="meshers", cells=(16, 8, 48))
    assert mesh.extract_surface(algorithm=None).n_open_edges == 0
    assert np.linalg.norm(mesh.points[:, :2], axis=1).min() >= 0.199
    assert mesh.volume == pytest.approx(np.pi * (0.4**2 - 0.2**2), rel=0.01)


def test_spherical_poles_are_not_silently_excised():
    shape = SphericalTpms(
        radius=2,
        surface_function=gyroid,
        offset=1.1,
        repeat_cell=(1, 0, 0),
        resolution=8,
    )
    with pytest.raises(ValueError, match="poles"):
        shape.generate_volume_mesh(backend="meshers")


def test_sector_periodicity_uses_world_rigid_transform(tmp_path):
    pytest.importorskip("meshers")
    shape = CylindricalTpms(
        radius=2,
        surface_function=gyroid,
        offset=1.1,
        resolution=8,
        center=(2, 3, 4),
        orientation=(13, 21, 7),
    )
    mesh = shape.generate_volume_mesh(
        backend="meshers", cells=16, periodic=(False, True, False)
    )
    pairs = mesh.field_data["periodic_pairs_y"]
    transform = mesh.field_data["periodic_transform_y"]
    expected = mesh.points[pairs[:, 0]] @ transform[:3, :3].T + transform[:3, 3]
    np.testing.assert_allclose(expected, mesh.points[pairs[:, 1]], atol=1e-10)
    assert "periodic_shift_y" not in mesh.field_data
    path = tmp_path / "sector.vtu"
    mesh.save(path)
    restored = pv.read(path)
    np.testing.assert_array_equal(
        restored.field_data["periodic_transform_y"], transform
    )


def test_bent_sweep_queries_reproduce_the_forward_map():
    theta = np.linspace(0, 0.6, 15)
    curve = np.column_stack(
        (3 * np.cos(theta), 3 * np.sin(theta), np.zeros_like(theta))
    )
    shape = Sweep(
        curve_points=curve,
        surface_function=gyroid,
        radial_max=0.2,
        radial_min=0.1,
        offset=1.1,
        resolution=8,
    )
    np.testing.assert_allclose(
        shape.raw_field(*shape.grid.points.T), shape.grid["surface"], atol=1e-9
    )


def test_grading_as_field_requires_a_fixed_reference():
    shape = Tpms(gyroid, offset=1.1, resolution=8)
    grading = NormedDistance(
        pv.Plane(direction=(1, 0, 0)), boundary_offset=1.2, furthest_offset=0.8
    )
    field = grading.as_field(shape.grid)
    points = shape.grid.points[[5, 17, 25]]
    np.testing.assert_allclose(field(*points[:1].T), field(*points.T)[:1])
