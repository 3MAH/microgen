import subprocess
import warnings

import numpy as np
import pytest
import pyvista as pv

from microgen import Tpms, is_periodic
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
    warnings.warn("MMG will not be used in these tests")
    USE_MMG = False


@pytest.fixture(name="box_mesh", scope="function")
def fixture_box_mesh() -> pv.UnstructuredGrid:
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
            [0.5, 0.0, 0.0],
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
    mesh = pv.UnstructuredGrid(elements, cell_types, nodes)

    return mesh


@pytest.fixture(name="gyroid_mesh", scope="function")
def fixture_gyroid_mesh() -> pv.UnstructuredGrid:
    gyroid_vtk = pv.UnstructuredGrid(
        Tpms(surface_function=gyroid, offset=1.0).generateVtk(type_part="sheet")
    )
    return gyroid_vtk


@pytest.fixture(name="non_periodic_mesh", scope="function")
def fixture_non_periodic_mesh() -> pv.UnstructuredGrid:
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
    mesh = pv.UnstructuredGrid(elements, cell_types, nodes)

    return mesh


@pytest.mark.parametrize(
    "shape",
    [
        "box_mesh",
        "gyroid_mesh",
    ],
)
def test_given_periodic_mesh_remesh_keeping_periodicity_for_fem_must_maintain_periodicity(
    shape: pv.UnstructuredGrid,
) -> None:
    if USE_MMG:
        # Act
        remeshed_shape = remesh_keeping_periodicity_for_fem(
            shape, hgrad=1.05
        ).to_pyvista()

        # Assert
        assert is_periodic(shape.points) and is_periodic(remeshed_shape.points)


def test_given_non_periodic_mesh_remesh_must_raise_inputmeshnotperiodicerror(
    non_periodic_mesh: pv.UnstructuredGrid,
) -> None:
    if USE_MMG:
        with pytest.raises(
            InputMeshNotPeriodicError, match="Input mesh is not periodic"
        ):
            remesh_keeping_periodicity_for_fem(non_periodic_mesh, hgrad=1.05)
