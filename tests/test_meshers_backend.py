"""Integration checks against the published meshers 0.1.0 wheel."""

import importlib.metadata

import numpy as np
import pytest
import pyvista as pv
from scipy.spatial.transform import Rotation

from microgen import Phase, Sphere, Tpms
from microgen.shape.implicit_ops import from_field
from microgen.shape.surface_functions import gyroid


def assert_tetrahedra(mesh):
    assert mesh.n_cells > 0
    assert np.all(mesh.celltypes == pv.CellType.TETRA)
    points = mesh.points[mesh.cells_dict[pv.CellType.TETRA]]
    determinants = np.linalg.det(points[:, 1:] - points[:, :1])
    assert np.all(determinants > 0)


def test_release_and_sphere_volume():
    assert importlib.metadata.version("meshers") == "0.1.0"
    mesh = Sphere(radius=1, center=(2, 3, 4)).generate_volume_mesh(resolution=16)
    assert_tetrahedra(mesh)
    assert mesh.volume == pytest.approx(4 * np.pi / 3, rel=0.04)
    assert mesh.center == pytest.approx((2, 3, 4), abs=0.02)


def test_surface_and_volume_share_boundary_and_do_not_alias():
    shape = from_field(
        lambda x, y, z: x * x + y * y + z * z - 1, bounds=(-1.2, 1.2) * 3
    )
    surface = shape.generate_surface_mesh(resolution=14)
    mesh = shape.generate_volume_mesh(resolution=14)
    assert surface.n_open_edges == 0
    assert surface.volume == pytest.approx(mesh.volume, rel=1e-10)
    surface.points[:] = 0
    assert shape.generate_surface_mesh(resolution=14).volume > 4


def test_phase_iso_and_exact_box_caps():
    phase = Phase(field=lambda x, y, z: x, bounds=(-1, 1) * 3, iso=0.25)
    mesh = phase.volume_mesh(resolution=10)
    assert_tetrahedra(mesh)
    assert mesh.volume == pytest.approx(5, abs=1e-10)
    boundary = phase.surface_mesh(resolution=10)
    assert boundary.n_open_edges == 0
    assert boundary.bounds == pytest.approx((-1, 0.25, -1, 1, -1, 1))


@pytest.mark.parametrize("part", ["sheet", "upper skeletal", "lower skeletal"])
def test_tpms_periodic_nodes_and_rigid_transform(part):
    rotation = Rotation.from_euler("z", 31, degrees=True)
    shape = Tpms(
        gyroid, offset=0.5, resolution=12, center=(2, 3, 4), orientation=rotation
    )
    mesh = shape.generate_meshers(part, periodic=(True,) * 3)
    for axis, pairs in enumerate(mesh.periodic_pairs):
        assert len(pairs)
        shift = rotation.apply(np.eye(3)[axis])
        np.testing.assert_allclose(
            mesh.points[pairs[:, 1]] - mesh.points[pairs[:, 0]],
            np.broadcast_to(shift, (len(pairs), 3)),
            atol=1e-10,
        )


def test_tpms_partition_and_offset_units():
    shape = Tpms(gyroid, offset=0.5, resolution=14)
    volumes = [
        shape.generate_volume_mesh(p).volume
        for p in ("sheet", "upper skeletal", "lower skeletal")
    ]
    assert 0.14 < volumes[0] < 0.19
    assert sum(volumes) == pytest.approx(1, abs=0.015)
    shape.offset = 1.0
    assert shape.generate_volume_mesh().volume > volumes[0] * 1.8


def test_callable_sheet_and_legacy_nodal_offset():
    shape = Tpms(gyroid, offset=lambda x, y, z: 0.6 + 0.1 * x, resolution=10)
    volume = shape.generate_volume_mesh().volume
    assert 0.15 < volume < 0.25
    shape.offset = np.asarray(shape.offset).copy()
    assert shape.generate_volume_mesh().volume == pytest.approx(volume, rel=0.04)


def test_empty_and_full_fields():
    for value, expected in [(1, 0), (-1, 1)]:
        shape = from_field(lambda x, y, z: value, bounds=(0, 1) * 3)
        mesh = shape.generate_volume_mesh(resolution=8)
        assert mesh.volume == pytest.approx(expected)


def test_invalid_resolution_and_mesher_failure_propagate():
    shape = Sphere()
    with pytest.raises(ValueError, match="resolution"):
        shape.generate_volume_mesh(resolution=1)
    with pytest.raises(ValueError):
        shape.generate_volume_mesh(resolution=10, minimum_quality=-1)


@pytest.mark.parametrize("name", ["CylindricalTpms", "SphericalTpms"])
def test_curved_tpms(name):
    import microgen

    shape = getattr(microgen, name)(
        radius=2,
        surface_function=gyroid,
        offset=1,
        resolution=10,
        repeat_cell=(1, 2, 2),
    )
    mesh = shape.generate_volume_mesh(optimize_passes=0)
    assert_tetrahedra(mesh)
    assert mesh.extract_surface(algorithm=None).n_open_edges == 0


def test_lattice_direct_periodic_meshing():
    from microgen import Cubic

    shape = Cubic(strut_radius=0.2)
    mesh = shape.generate_volume_mesh(resolution=10, optimize_passes=0)
    assert_tetrahedra(mesh)
    for axis in "xyz":
        assert len(mesh.field_data[f"periodic_pairs_{axis}"])


def test_density_target():
    shape = Tpms(gyroid, density=0.3, resolution=10)
    mesh = shape.generate_volume_mesh()
    assert mesh.volume == pytest.approx(0.3, abs=0.003)


def test_full_angular_chart_retains_legacy_path():
    from microgen import CylindricalTpms

    shape = CylindricalTpms(
        radius=1,
        surface_function=gyroid,
        offset=0.5,
        resolution=8,
        repeat_cell=(1, 0, 1),
    )
    assert shape.generate_volume_mesh().n_cells > 0
    with pytest.raises(NotImplementedError, match="parametric mesher"):
        shape.generate_meshers()


def test_infill_stays_inside_envelope():
    from microgen import Infill

    envelope = pv.Sphere(radius=0.8, theta_resolution=12, phi_resolution=12)
    shape = Infill(
        obj=envelope, surface_function=gyroid, offset=0.8, repeat_cell=1, resolution=10
    )
    mesh = shape.generate_volume_mesh(optimize_passes=0)
    assert_tetrahedra(mesh)
    assert 0 < mesh.volume < envelope.volume
    assert mesh.extract_surface(algorithm=None).n_open_edges == 0
