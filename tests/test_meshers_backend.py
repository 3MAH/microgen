"""Direct TPMS volume meshing with the published meshers wheel."""

import importlib.metadata
import json
import meshers
import numpy as np
import pytest
import pyvista as pv
from scipy.spatial.transform import Rotation
from microgen import CylindricalTpms, Tpms
from microgen.shape.surface_functions import gyroid


@pytest.mark.parametrize("part", ["sheet", "upper skeletal", "lower skeletal"])
def test_periodic_quality_and_boundary_triangulations(part):
    assert importlib.metadata.version("meshers") == "0.1.0"
    rotation = Rotation.from_euler("z", 31, degrees=True)
    shape = Tpms(
        gyroid, offset=0.5, resolution=16, center=(2, 3, 4), orientation=rotation
    )
    native = shape.generate_meshers(part, periodic=(True,) * 3)
    assert native.diagnostics["minimum_mmg_quality"] >= 0.1
    assert native.diagnostics["sampled_surface_error"] <= 0.01
    vertices = native.points[native.tetrahedra]
    assert np.all(np.linalg.det(vertices[:, 1:] - vertices[:, :1]) > 0)
    for axis, pairs in enumerate(native.periodic_pairs):
        assert len(pairs)
        shift = rotation.apply(np.eye(3)[axis])
        np.testing.assert_allclose(
            native.points[pairs[:, 1]] - native.points[pairs[:, 0]],
            np.broadcast_to(shift, (len(pairs), 3)),
            atol=1e-10,
        )
        lookup = dict(pairs.tolist())
        low = native.surface[native.boundary_tags == 2 * axis + 1]
        high = native.surface[native.boundary_tags == 2 * axis + 2]
        assert {tuple(sorted(lookup[int(i)] for i in f)) for f in low} == {
            tuple(sorted(f)) for f in high
        }


def test_volume_quality_metadata_and_legacy_surface_preserved():
    shape = Tpms(gyroid, offset=0.5, resolution=16)
    old_points = shape.generate_surface_mesh().points.copy()
    grid = shape.generate_volume_mesh(periodic=(True,) * 3)
    assert np.all(grid.celltypes == pv.CellType.TETRA)
    assert grid.cell_data["MMGQuality"].min() >= 0.1
    assert grid.cell_data["Volume"].min() > 0
    assert grid.volume == pytest.approx(grid.cell_data["Volume"].sum())
    assert (
        json.loads(grid.field_data["meshers_diagnostics"][0])["minimum_mmg_quality"]
        >= 0.1
    )
    assert grid.extract_surface(algorithm=None).n_open_edges == 0
    np.testing.assert_array_equal(shape.generate_surface_mesh().points, old_points)


def test_rejection_does_not_fall_back_to_legacy(monkeypatch):
    shape = Tpms(gyroid, offset=0.5, resolution=12)

    def forbidden(*args, **kwargs):
        pytest.fail("quality failure must not invoke the legacy mesher")

    monkeypatch.setattr(shape, "_generate_legacy_volume_mesh", forbidden)
    with pytest.raises(meshers.MeshingError, match="quality"):
        shape.generate_volume_mesh(minimum_quality=0.999)


def test_nonperiodic_grading_rejected_as_periodic():
    shape = Tpms(gyroid, offset=lambda x, y, z: 0.6 + 0.2 * x, resolution=12)
    with pytest.raises(meshers.MeshingError, match="periodic"):
        shape.generate_volume_mesh(periodic=(True,) * 3)


def test_native_density_calibration():
    shape = Tpms(gyroid, density=0.3, resolution=16)
    grid = shape.generate_volume_mesh(periodic=(True,) * 3)
    assert grid.volume == pytest.approx(0.3, abs=0.003)
    assert shape.density == 0.3


def test_density_failure_restores_input():
    shape = Tpms(gyroid, density=0.3, resolution=10)
    before = shape.offset
    with pytest.raises(meshers.MeshingError):
        shape.generate_meshers(minimum_quality=0.999)
    assert shape.density == 0.3
    assert shape.offset == before


def test_retained_full_wrap_is_explicit():
    shape = CylindricalTpms(
        radius=1,
        surface_function=gyroid,
        offset=0.5,
        resolution=8,
        repeat_cell=(1, 0, 1),
    )
    with pytest.warns(RuntimeWarning, match="legacy clipped grid"):
        assert shape.generate_volume_mesh().n_cells > 0
    with pytest.raises(NotImplementedError):
        shape.generate_meshers()
    with pytest.raises(NotImplementedError):
        shape.generate_volume_mesh(minimum_quality=0.1)


def test_resolution_limit_includes_repeats():
    shape = Tpms(gyroid, offset=0.5, resolution=10, repeat_cell=(14, 1, 1))
    with pytest.raises(ValueError, match="129"):
        shape.generate_volume_mesh()


def test_backend_selector_is_not_an_option():
    shape = Tpms(gyroid, offset=0.5, resolution=10)
    with pytest.raises(TypeError, match="backend"):
        shape.generate_volume_mesh(backend="mmg")
