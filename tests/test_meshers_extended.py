"""Field and domain controls for sampled, mapped and clipped TPMS."""

import numpy as np
import pytest
import pyvista as pv

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
    corners = mesh.points[mesh.cells_dict[np.uint8(pv.CellType.TETRA)]]
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


@pytest.mark.parametrize("graded", [False, True])
def test_adaptive_sweep_keeps_both_sheet_constraints(monkeypatch, graded):
    meshers = pytest.importorskip("meshers")
    # Keep the legacy Cartesian compiler contract covered with newer wheels.
    monkeypatch.setattr(meshers, "ADAPTIVE_MAPPED_INTERSECTIONS", False, raising=False)
    shape = Sweep(
        curve_points=np.array([[0.0, 0, 0], [0.0, 0, 1]]),
        surface_function=gyroid,
        radial_max=0.4,
        radial_min=0.2,
        offset=(lambda x, y, z: 1.1 + 0.1 * z) if graded else 1.1,
        resolution=8,
    )
    points = np.random.default_rng(15).uniform(
        np.asarray(shape.bounds)[::2], np.asarray(shape.bounds)[1::2], (128, 3)
    )

    def checked(field, **options):
        compiled = meshers.compile_field(field)
        # Retention is the native compiler contract used by adaptive meshing.
        branches = compiled._intersection_fields()
        assert branches is not None and len(branches) == 6
        half = (1.1 + 0.1 * points[:, 2] if graded else 1.1) / 2
        expected = np.maximum(
            np.abs(shape.raw_field(*points.T)) - half,
            shape._cell_box().require_field()(*points.T),
        )
        np.testing.assert_allclose(compiled.values(points), expected, atol=1e-12)
        np.testing.assert_allclose(
            np.max([branch.values(points) for branch in branches], axis=0),
            expected,
            atol=1e-12,
        )
        raise RuntimeError("checked adaptive constraints")

    monkeypatch.setattr(meshers, "generate", checked)
    with pytest.raises(RuntimeError, match="checked adaptive constraints"):
        shape.generate_volume_mesh(backend="meshers", adapt=True)


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


@pytest.mark.parametrize(
    "cells, passes", [((16, 8, 48), 4), ((24, 12, 64), 4), (16, 1)]
)
def test_annular_sweep_sheet_adapts_with_named_physical_constraints(cells, passes):
    meshers = pytest.importorskip("meshers")
    if not getattr(meshers, "ADAPTIVE_MAPPED_INTERSECTIONS", False):
        pytest.skip("requires named seed adaptation support")
    shape = Sweep(
        curve_points=np.array([[0.0, 0, 0], [0.0, 0, 0.5], [0.0, 0, 1]]),
        surface_function=gyroid,
        radial_max=0.4,
        radial_min=0.2,
        offset=1.1,
        resolution=8,
    )
    mesh = shape.generate_volume_mesh(
        backend="meshers",
        cells=cells,
        adapt=True,
        adapt_tolerance=0.0018,
        geometry_tolerance=0.002,
        min_length=0.002,
        max_length=0.08,
        gradation=0.6,
        max_tetrahedra=250_000,
        optimize_passes=passes,
    )
    assert mesh.n_cells <= 250_000 and mesh.n_points <= 250_000
    assert mesh.extract_surface(algorithm=None).n_open_edges == 0
    assert mesh.volume == pytest.approx(0.13433935, rel=0.01)
    corners = mesh.points[mesh.cells_dict[np.uint8(pv.CellType.TETRA)]]
    det = np.linalg.det(corners[:, 1:] - corners[:, :1])
    edges = sum(
        np.sum((corners[:, a] - corners[:, b]) ** 2, axis=1)
        for a in range(4)
        for b in range(a + 1, 4)
    )
    assert np.all(det > 0) and np.sqrt(432 * det**2 / edges**3).min() >= 0.1
    assert len(mesh.field_data["meshers_constraint_names"]) == 6


@pytest.mark.parametrize("direction", [(0, 0, 1), (1, 2, 3), (-2, 0.3, 1), (1, 0, 0)])
def test_straight_sweep_map_matches_frames_and_compiles(direction):
    meshers = pytest.importorskip("meshers")
    curve = np.array([0.2, -0.4, 0.7]) + np.array([0, 0.13, 1])[:, None] * direction
    shape = Sweep(
        curve_points=curve,
        surface_function=gyroid,
        radial_min=0.2,
        radial_max=0.4,
        cell_size=(0.8, 1.2, 0.9),
        repeat_cell=(2, 1, 3),
        offset=1.1,
        resolution=8,
    )
    extent = shape.cell_size * shape.repeat_cell
    query = (
        np.random.default_rng(423).uniform(-0.7, 0.7, (3, 3, 7)) * extent[:, None, None]
    )
    s = np.clip(
        (query[0] / extent[0] + 0.5) * shape._curve_total_length,
        0,
        shape._curve_total_length,
    )
    r = shape.radial_min + (query[1] / extent[1] + 0.5) * (
        shape.radial_max - shape.radial_min
    )
    theta = (query[2] / extent[2] + 0.5) * 2 * np.pi
    c, normal, tangent = shape._curve_frame_at(s.ravel(order="F"))
    expected = c + r.ravel(order="F")[:, None] * (
        np.cos(theta.ravel(order="F"))[:, None] * normal
        + np.sin(theta.ravel(order="F"))[:, None] * np.cross(tangent, normal)
    )
    actual = np.column_stack(
        [a.ravel(order="F") for a in shape._map_coordinates(*query)]
    )
    np.testing.assert_allclose(actual, expected, atol=2e-13, rtol=5e-14)
    points = np.column_stack([a.ravel(order="F") for a in query])
    for component in range(3):

        def coordinate(x, y, z, component=component):
            return shape._map_coordinates(x, y, z)[component]

        compiled = meshers.compile_field(coordinate)
        np.testing.assert_allclose(
            compiled.values(points), expected[:, component], atol=2e-13, rtol=5e-14
        )


def test_curved_sweep_envelope_retains_interior_and_local_end_caps():
    theta = np.linspace(0, 1.5 * np.pi, 61)
    curve = np.column_stack(
        (3 * np.cos(theta), 3 * np.sin(theta), np.zeros_like(theta))
    )
    shape = Sweep(
        curve, gyroid, radial_min=0.1, radial_max=0.2, offset=1.1, resolution=4
    )
    parameter = np.column_stack(
        (np.linspace(-0.45, 0.45, 13), np.zeros(13), np.zeros(13))
    )
    points = np.column_stack(shape._map_coordinates(*parameter.T))
    envelope = shape._cell_box().require_field()
    assert np.all(envelope(*points.T) < 0)
    centers, normals, tangents = shape._curve_frame_at(
        np.array([0.0, shape._curve_total_length])
    )
    caps = centers + 0.15 * normals
    outside = caps + np.array([-0.05, 0.05])[:, None] * tangents
    inside = caps - np.array([-0.05, 0.05])[:, None] * tangents
    assert np.all(envelope(*outside.T) > 0)
    assert np.all(envelope(*inside.T) < 0)
    radial_outside = centers + 0.25 * normals
    assert np.all(envelope(*radial_outside.T) > 0)


@pytest.mark.parametrize("part", ["sheet", "lower skeletal", "upper skeletal"])
def test_full_density_sweep_adaptive_fallback_uses_only_envelope(monkeypatch, part):
    meshers = pytest.importorskip("meshers")
    shape = Sweep(
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
        gyroid,
        radial_min=0.1,
        radial_max=0.2,
        density=1,
        resolution=4,
    )
    parameter = np.array([[0.125, 0.125, 0.125]])
    points = np.column_stack(shape._map_coordinates(*parameter.T))

    def checked(field, **options):
        assert options["adapt"] and options["order"] == 2
        envelope = shape._cell_box().require_field()
        expected = envelope(*points.T)
        assert np.all(expected < 0)
        np.testing.assert_array_equal(field(*points.T), expected)
        raise RuntimeError("checked full-density fallback")

    monkeypatch.setattr(meshers, "generate", checked)
    with pytest.raises(RuntimeError, match="checked full-density fallback"):
        shape.generate_volume_mesh(
            type_part=part, backend="meshers", adapt=True, order=2, cells=4
        )
