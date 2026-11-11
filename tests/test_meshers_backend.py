"""Opt-in volume meshing preserves fields, placement and periodic connectivity."""

import json
import sys

import numpy as np
import pytest
import pyvista as pv
from scipy.spatial.transform import Rotation

from microgen import CylindricalTpms, Sphere, Spinodoid, Tpms
from microgen.shape.implicit_ops import from_field
from microgen.shape.surface_functions import gyroid, schwarz_p
from microgen.shape.tpms_grading import OffsetGrading


@pytest.fixture
def native_backend():
    return pytest.importorskip("meshers")


def assert_native_grid(grid, order=1):
    assert isinstance(grid, pv.UnstructuredGrid) and grid.n_cells > 0
    cell_type = pv.CellType.TETRA if order == 1 else pv.CellType.QUADRATIC_TETRA
    assert np.all(grid.celltypes == cell_type)
    assert grid.cell_data["MMGQuality"].min() >= 0.1
    assert grid.cell_data["CornerSignedVolume"].min() > 0
    assert grid.extract_surface(algorithm=None).n_open_edges == 0
    diagnostics = json.loads(grid.field_data["meshers_diagnostics"][0])
    assert diagnostics["evaluator"] == "compiled"
    assert diagnostics["threads"] == 1


def assert_periodic_faces(grid, rotation, lengths):
    faces = grid.field_data["meshers_boundary_faces"]
    tags = grid.field_data["meshers_boundary_tags"]
    for axis, name in enumerate("xyz"):
        pairs = grid.field_data[f"periodic_pairs_{name}"]
        expected = rotation.apply(np.eye(3)[axis] * lengths[axis])
        np.testing.assert_allclose(
            grid.points[pairs[:, 1]] - grid.points[pairs[:, 0]],
            np.broadcast_to(expected, (len(pairs), 3)),
            atol=1e-11,
        )
        np.testing.assert_allclose(
            grid.field_data[f"periodic_shift_{name}"][0], expected
        )
        lookup = dict(pairs.tolist())
        low, high = faces[tags == 2 * axis + 1], faces[tags == 2 * axis + 2]
        assert len(low) > 0 and len(low) == len(high)
        assert {tuple(sorted(lookup[int(i)] for i in f)) for f in low} == {
            tuple(sorted(f)) for f in high
        }


@pytest.mark.parametrize("kind", ["shape", "tpms", "spinodoid"])
def test_vtk_default_does_not_require_meshers(monkeypatch, kind):
    shapes = {
        "shape": lambda: Sphere(radius=0.35),
        "tpms": lambda: Tpms(gyroid, offset=1.1, resolution=8),
        "spinodoid": lambda: Spinodoid(offset=0, resolution=8, seed=4),
    }
    shape = shapes[kind]()
    monkeypatch.setitem(sys.modules, "meshers", None)
    original = shape.generate_volume_mesh(resolution=8)
    selected = shape.generate_volume_mesh(backend="vtk", resolution=8)
    np.testing.assert_array_equal(original.points, selected.points)
    np.testing.assert_array_equal(original.cells, selected.cells)
    with pytest.raises(ImportError, match="optional Meshers package"):
        shape.generate_volume_mesh(backend="meshers", cells=8)
    with pytest.raises(ValueError, match="Unknown volume-mesh backend"):
        shape.generate_volume_mesh(backend="unknown")


@pytest.mark.parametrize("order", [1, 2])
def test_primitive_world_field_is_not_translated_twice(native_backend, order):
    shape = Sphere(radius=0.35, center=(2, 3, 4))
    grid = shape.generate_volume_mesh(
        backend="meshers", cells=12, order=order, optimize_passes=1
    )
    assert_native_grid(grid, order)
    np.testing.assert_allclose(grid.center, shape.center, atol=1e-4)
    assert grid.cells_dict[
        pv.CellType.TETRA if order == 1 else pv.CellType.QUADRATIC_TETRA
    ].shape[1] == (4 if order == 1 else 10)
    assert "Volume" not in grid.cell_data  # Corner volume is labeled explicitly.


@pytest.mark.parametrize("order", [1, 2])
@pytest.mark.parametrize("part", ["sheet", "lower skeletal", "upper skeletal"])
def test_tpms_part_units_placement_and_periodic_facets(native_backend, part, order):
    rotation = Rotation.from_euler("z", 31, degrees=True)
    shape = Tpms(
        gyroid if part == "sheet" else schwarz_p,
        offset=1.1 if part == "sheet" else 0.0,
        resolution=12,
        center=(2, 3, 4),
        orientation=rotation,
    )
    original = shape.generate_volume_mesh()
    grid = shape.generate_volume_mesh(
        type_part=part,
        backend="meshers",
        cells=12,
        order=order,
        periodic=(True,) * 3,
        geometry_tolerance=0.01,
        optimize_passes=1,
        max_tetrahedra=250_000,
    )
    assert_native_grid(grid, order)
    assert_periodic_faces(grid, rotation, (1, 1, 1))
    if part != "sheet":
        assert grid.volume == pytest.approx(0.5, abs=0.015)
    unchanged = shape.generate_volume_mesh()
    np.testing.assert_array_equal(original.points, unchanged.points)
    np.testing.assert_array_equal(original.cells, unchanged.cells)


def test_spinodoid_snapshot_matches_existing_field_and_freezes_coefficients():
    from microgen._meshers import spinodoid_field

    shape = Spinodoid(
        k0=4 * np.pi,
        bandwidth=3,
        offset=0.05,
        resolution=16,
        cell_size=(1, 1.25, 0.75),
        repeat_cell=(2, 1, 1),
        seed=20261026,
    )
    points = np.random.default_rng(4).uniform(0, 1, (4096, 3))
    points *= shape.cell_size * shape.repeat_cell
    field = spinodoid_field(shape._frep)
    expected = shape.evaluate(*points.T)
    np.testing.assert_allclose(field(*points.T), expected, atol=1e-12, rtol=1e-12)
    np.testing.assert_array_equal(field(*points.T) < 0, expected < 0)
    shape._frep.amplitudes *= 2
    np.testing.assert_allclose(field(*points.T), expected, atol=1e-12, rtol=1e-12)


def test_spinodoid_meshing_keeps_fft_realization_and_local_origin(native_backend):
    rotation = Rotation.from_euler("z", 31, degrees=True)
    shape = Spinodoid(
        k0=4 * np.pi,
        bandwidth=3,
        offset=0,
        resolution=16,
        anisotropy=(1, 1, 2),
        seed=20261026,
        mode_amplitude_threshold=0.05,
        center=(2, 3, 4),
        orientation=rotation,
    )
    coefficients = np.column_stack(
        (shape._frep.modes, shape._frep.amplitudes, shape._frep.phases)
    ).copy()
    assert "grid" not in shape.__dict__
    grid = shape.generate_volume_mesh(
        backend="meshers",
        cells=16,
        periodic=(True,) * 3,
        adapt=True,
        adapt_tolerance=0.0018,
        geometry_tolerance=0.002,
        min_length=0.002,
        max_length=0.1,
        gradation=0.6,
        optimize_passes=1,
        max_tetrahedra=250_000,
    )
    assert_native_grid(grid)
    assert_periodic_faces(grid, rotation, (1, 1, 1))
    local = rotation.inv().apply(grid.points - np.asarray(shape.center))
    np.testing.assert_allclose(local.min(axis=0), 0, atol=1e-11)
    np.testing.assert_allclose(local.max(axis=0), 1, atol=1e-11)
    assert shape.resolution == 16 and "grid" not in shape.__dict__
    np.testing.assert_array_equal(
        np.column_stack(
            (shape._frep.modes, shape._frep.amplitudes, shape._frep.phases)
        ),
        coefficients,
    )


def test_composed_field_and_cells_override_resolution(native_backend, monkeypatch):
    shape = Sphere(radius=0.3, center=(-0.1, 0, 0)) | Sphere(
        radius=0.3, center=(0.1, 0, 0)
    )

    def no_grid(*args, **kwargs):
        pytest.fail("Meshers must use the exact callable, not a sampled VTK grid")

    monkeypatch.setattr(shape, "_sample_implicit_grid", no_grid)
    grid = shape.generate_volume_mesh(
        backend="meshers", resolution=200, cells=16, optimize_passes=1
    )
    assert_native_grid(grid)


def test_native_quality_error_propagates_without_vtk_retry(native_backend, monkeypatch):
    shape = Sphere(radius=0.35)
    calls = []
    generate = native_backend.generate

    def record(*args, **kwargs):
        calls.append(kwargs)
        return generate(*args, **kwargs)

    monkeypatch.setattr(native_backend, "generate", record)
    monkeypatch.setattr(
        shape, "_sample_implicit_grid", lambda *a, **k: pytest.fail("fallback")
    )
    with pytest.raises(native_backend.MeshingError, match="quality"):
        shape.generate_volume_mesh(
            backend="meshers", cells=8, geometry_tolerance=1, minimum_quality=0.999
        )
    assert len(calls) == 1


def test_nonperiodic_grading_is_not_canonicalized(native_backend):
    shape = Tpms(gyroid, offset=lambda x, y, z: 1.1 + 0.2 * z, resolution=12)
    with pytest.raises(ValueError, match="periodic"):
        shape.generate_volume_mesh(backend="meshers", cells=12, periodic=(True,) * 3)


def test_existing_density_calibration_is_reused(native_backend):
    shape = Tpms(gyroid, density=0.3, resolution=16)
    shape.generate_volume_mesh()
    offset = shape.offset
    grid = shape.generate_volume_mesh(backend="meshers", cells=16, optimize_passes=1)
    assert_native_grid(grid)
    assert grid.volume == pytest.approx(0.3, abs=0.015)
    assert shape.offset == offset and shape.density == 0.3


def test_repeated_full_cell_and_saved_periodic_metadata(native_backend, tmp_path):
    rotation = Rotation.from_euler("z", 31, degrees=True)
    shape = Tpms(
        gyroid,
        density=1,
        cell_size=(1, 1.25, 0.75),
        repeat_cell=(2, 1, 1),
        resolution=8,
        center=(2, 3, 4),
        orientation=rotation,
    )
    grid = shape.generate_volume_mesh(
        backend="meshers", order=2, periodic=(True,) * 3, optimize_passes=1
    )
    assert_native_grid(grid, order=2)
    assert grid.volume == pytest.approx(2 * 1.25 * 0.75, abs=1e-10)
    assert_periodic_faces(grid, rotation, (2, 1.25, 0.75))
    path = tmp_path / "periodic.vtu"
    grid.save(path)
    restored = pv.read(path)
    assert_periodic_faces(restored, rotation, (2, 1.25, 0.75))
    for name in grid.field_data:
        np.testing.assert_array_equal(restored.field_data[name], grid.field_data[name])


def test_sampled_offsets_require_an_explicit_policy():
    class SampledOffset(OffsetGrading):
        def compute_offset(self, grid):
            return np.full(grid.n_points, 1.1)

    shape = Tpms(gyroid, offset=SampledOffset(), resolution=8)
    assert shape.generate_volume_mesh().n_cells
    with pytest.raises(NotImplementedError, match="interpolation policy"):
        shape.generate_volume_mesh(backend="meshers", cells=8)


def test_mapped_tpms_are_explicitly_unsupported():
    shape = CylindricalTpms(radius=2, surface_function=gyroid, offset=1.1, resolution=8)
    with pytest.raises(NotImplementedError, match="domain adapter"):
        shape.generate_volume_mesh(backend="meshers", cells=8)


def test_shape_region_cannot_be_changed_through_backend_options(native_backend):
    shape = Sphere(radius=0.35)
    with pytest.raises(ValueError, match="solid region"):
        shape.generate_volume_mesh(backend="meshers", cells=8, band=(-0.1, 0.1))


def test_empty_field_returns_an_empty_grid(native_backend):
    shape = from_field(lambda x, y, z: x * 0 + 1, bounds=(-1, 1, -1, 1, -1, 1))
    grid = shape.generate_volume_mesh(backend="meshers", cells=4)
    assert isinstance(grid, pv.UnstructuredGrid) and grid.n_cells == grid.n_points == 0
