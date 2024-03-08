import subprocess

import numpy as np
import pytest
import pyvista as pv

from microgen import BoxMesh, Tpms, is_periodic
from microgen.remesh import (
    InputMeshNotPeriodicError,
    remesh_keeping_periodicity_for_fem,
)
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

    return mesh


@pytest.fixture(name="gyroid_mesh", scope="function")
def fixture_gyroid_mesh() -> BoxMesh:
    gyroid_vtk = pv.UnstructuredGrid(
        Tpms(surface_function=gyroid, offset=1.0).generateVtk(type_part="sheet")
    )
    gyroid_mesh = BoxMesh.from_pyvista(gyroid_vtk)
    return gyroid_mesh


@pytest.fixture(name="non_periodic_box_mesh", scope="function")
def fixture_non_periodic_box_mesh() -> BoxMesh:
    nodes = np.array(
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
            [0.5, 0.1, 0.0],
            [-0.5, 0.0, 0.0],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, -0.5],
        ]
    )
    elements = np.array(
        [
            [4, 11, 6, 0, 14],
            [4, 0, 9, 1, 11],
            [4, 3, 10, 0, 13],
            [4, 9, 0, 4, 11],
            [4, 5, 11, 0, 13],
            [4, 0, 11, 1, 13],
            [4, 0, 10, 5, 13],
            [4, 11, 5, 1, 13],
            [4, 3, 10, 2, 12],
            [4, 6, 10, 0, 14],
            [4, 1, 9, 4, 11],
            [4, 10, 3, 0, 12],
            [4, 10, 3, 5, 13],
            [4, 0, 10, 2, 14],
            [4, 2, 10, 0, 12],
            [4, 4, 11, 0, 14],
            [4, 6, 11, 4, 14],
            [4, 10, 6, 2, 14],
            [4, 5, 10, 0, 11],
            [4, 0, 10, 6, 11],
            [4, 10, 5, 6, 11],
            [4, 7, 9, 0, 12],
            [4, 0, 9, 8, 12],
            [4, 9, 7, 8, 12],
            [4, 7, 9, 1, 13],
            [4, 1, 9, 0, 13],
            [4, 9, 7, 0, 13],
            [4, 0, 12, 3, 13],
            [4, 3, 12, 7, 13],
            [4, 12, 0, 7, 13],
            [4, 8, 12, 2, 14],
            [4, 2, 12, 0, 14],
            [4, 12, 8, 0, 14],
            [4, 8, 9, 0, 14],
            [4, 0, 9, 4, 14],
            [4, 9, 8, 4, 14],
        ]
    )
    cell_types = np.full(elements.shape[0], pv.CellType.TETRA, dtype=np.uint8)
    grid = pv.UnstructuredGrid(elements, cell_types, nodes)
    mesh = BoxMesh.from_pyvista(grid)

    return mesh


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
) -> None:
    if USE_MMG:
        # Arrange

        vtk_mesh = request.getfixturevalue(shape).to_pyvista()

        # Act

        remeshed_vtk_mesh = remesh_keeping_periodicity_for_fem(shape, hgrad=1.05)

        # Assert
        assert is_periodic(vtk_mesh.points) and is_periodic(remeshed_vtk_mesh.points)


def test_given_non_periodic_box_mesh_remesh_must_raise_inputmeshnotperiodicerror(
    non_periodic_box_mesh: BoxMesh,
) -> None:
    if USE_MMG:
        with pytest.raises(
            InputMeshNotPeriodicError, match="Input mesh is not periodic"
        ):
            remesh_keeping_periodicity_for_fem(non_periodic_box_mesh, hgrad=1.05)
