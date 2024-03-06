import subprocess
from pathlib import Path

import gmsh
import numpy as np
import numpy.typing as npt
import pytest
import pyvista as pv
from microgen import BoxMesh, Rve, Tpms, is_periodic, remesh
from microgen.shape.surface_functions import gyroid

_MESH_DIM = 3
_BOUNDARY_DIM = 2

USE_MMG = False
try:
    subprocess.check_output("mmg3d_O3", stderr=subprocess.STDOUT)
    USE_MMG = True
except (subprocess.CalledProcessError, FileNotFoundError):
    print("mmg command did not work, check if it is installed or contact a developer")
    USE_MMG = False


def _get_mesh_nodes_coords(mesh_name: str) -> npt.NDArray[np.float_]:
    gmsh.initialize()
    gmsh.open(mesh_name)

    _, nodes_coords, _ = gmsh.model.mesh.getNodes()
    gmsh.finalize()

    nodes_coords = nodes_coords.reshape(-1, _MESH_DIM)

    return nodes_coords


@pytest.fixture(name="tmp_dir", scope="function")
def fixture_tmp_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    tmp_dir_name = "test_tmp_dir"
    return tmp_path_factory.mktemp(tmp_dir_name)


@pytest.fixture(name="tmp_mesh_filename", scope="function")
def fixture_tmp_mesh_filename(tmp_dir: Path) -> str:
    return (tmp_dir / "shape.mesh").as_posix()


@pytest.fixture(name="tmp_output_mesh_filename", scope="function")
def fixture_tmp_output_mesh_filename(tmp_dir: Path) -> str:
    return (tmp_dir / "shape.o.mesh").as_posix()


def _default_rve() -> Rve:
    return Rve(dim_x=1.0, dim_y=1.0, dim_z=1.0, center=(0.5, 0.5, 0.5))


@pytest.fixture(name="box_mesh", scope="function")
def fixture_box_mesh() -> BoxMesh:
    nodes_array = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [-0.5, -0.5, -0.5],
            [-0.5, -0.5, 0.5],
            [0.5, 0.5, -0.5],
            [0.5, -0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, -0.5],
            [0.0, 0.5, 0.0],
            [0.0, -0.5, 0.0],
            [0.5, 0.0, 0.0],
            [-0.5, 0.0, 0.0],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, -0.5],
        ]
    )

    elements_dict = {
        pv.CellType.TETRA: np.array(
            [
                [11, 6, 0, 14],
                [0, 9, 1, 11],
                [3, 10, 0, 13],
                [9, 0, 4, 11],
                [5, 11, 0, 13],
                [0, 11, 1, 13],
                [0, 10, 5, 13],
                [11, 5, 1, 13],
                [3, 10, 2, 12],
                [6, 10, 0, 14],
                [1, 9, 4, 11],
                [10, 3, 0, 12],
                [10, 3, 5, 13],
                [0, 10, 2, 14],
                [2, 10, 0, 12],
                [4, 11, 0, 14],
                [6, 11, 4, 14],
                [10, 6, 2, 14],
                [5, 10, 0, 11],
                [0, 10, 6, 11],
                [10, 5, 6, 11],
                [7, 9, 0, 12],
                [0, 9, 8, 12],
                [9, 7, 8, 12],
                [7, 9, 1, 13],
                [1, 9, 0, 13],
                [9, 7, 0, 13],
                [0, 12, 3, 13],
                [3, 12, 7, 13],
                [12, 0, 7, 13],
                [8, 12, 2, 14],
                [2, 12, 0, 14],
                [12, 8, 0, 14],
                [8, 9, 0, 14],
                [0, 9, 4, 14],
                [9, 8, 4, 14],
            ]
        )
    }

    mesh = BoxMesh(nodes_array, elements_dict)
    mesh.rve = _default_rve()

    return mesh


@pytest.fixture(name="gyroid_mesh", scope="function")
def fixture_gyroid_mesh() -> BoxMesh:
    gyroid_vtk = pv.UnstructuredGrid(
        Tpms(surface_function=gyroid, offset=1.0).generateVtk(type_part="sheet")
    )
    gyroid_mesh = BoxMesh.from_pyvista(gyroid_vtk)
    return gyroid_mesh


@pytest.mark.parametrize(
    "shape",
    [
        "box_mesh",
        "gyroid_mesh",
    ],
)
def test_given_periodic_mesh_remesh_keeping_periodicity_for_fem_must_maintain_periodicity(
    shape: BoxMesh,
    request,
    tmp_mesh_filename: str,
    tmp_output_mesh_filename: str,
) -> None:
    if USE_MMG:
        # Arrange

        vtk_mesh = request.getfixturevalue(shape).to_pyvista()
        pv.save_meshio(tmp_mesh_filename, vtk_mesh)

        # Act

        remesh.remesh_keeping_periodicity_for_fem(
            shape, tmp_output_mesh_filename, hgrad=1.05
        )

        # Assert
        assert is_periodic(_get_mesh_nodes_coords(tmp_mesh_filename)) and is_periodic(
            _get_mesh_nodes_coords(tmp_output_mesh_filename)
        )
