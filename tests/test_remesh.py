"""Test for the remesh module."""

import shutil
import warnings

import numpy as np
import pytest
import pyvista as pv
from _pytest.fixtures import FixtureRequest

from microgen import BoxMesh, Tpms, is_periodic
from microgen.remesh import (
    InputMeshNotPeriodicError,
    remesh_keeping_periodicity_for_fem,
)
from microgen.shape.surface_functions import gyroid

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/
# ruff: noqa: E501 line-too-long https://docs.astral.sh/ruff/rules/line-too-long/


_MESH_DIM = 3
_BOUNDARY_DIM = 2

USE_MMG = shutil.which("mmg3d_O3") is not None
if not USE_MMG:
    warnings.warn("MMG will not be used in these tests", stacklevel=2)


@pytest.fixture(name="box_mesh")
def fixture_box_mesh() -> BoxMesh:
    """Return a periodic box mesh."""
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
        ],
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
            ],
        ),
    }

    return BoxMesh(nodes_array, elements_dict)


@pytest.fixture(name="gyroid_mesh")
def fixture_gyroid_mesh() -> pv.UnstructuredGrid:
    """Return a gyroid mesh."""
    return Tpms(surface_function=gyroid, offset=1.0).grid_sheet


@pytest.fixture(name="non_periodic_mesh")
def fixture_non_periodic_mesh() -> pv.UnstructuredGrid:
    """Return a non-periodic mesh."""
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
        ],
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
        ],
    )
    cell_types = np.full(elements.shape[0], pv.CellType.TETRA, dtype=np.uint8)
    return pv.UnstructuredGrid(elements, cell_types, nodes)


@pytest.mark.parametrize(
    "shape",
    [
        "box_mesh",
        "gyroid_mesh",
    ],
)
def test_given_periodic_mesh_remesh_keeping_periodicity_for_fem_must_maintain_periodicity(
    shape: str,
    request: FixtureRequest,
) -> None:
    """Test that remesh_keeping_periodicity_for_fem maintains periodicity."""
    # Arrange
    input_mesh = request.getfixturevalue(shape)
    # Act
    if USE_MMG:
        edge_length_gradient = 1.05
        remeshed_shape = remesh_keeping_periodicity_for_fem(
            input_mesh,
            hgrad=edge_length_gradient,
        )

        if isinstance(input_mesh, BoxMesh):
            input_mesh = input_mesh.to_pyvista()
            remeshed_shape = remeshed_shape.to_pyvista()
        # Assert
        assert is_periodic(input_mesh.points)
        assert is_periodic(remeshed_shape.points)


def test_given_non_periodic_mesh_remesh_must_raise_inputmeshnotperiodicerror(
    non_periodic_mesh: pv.UnstructuredGrid,
) -> None:
    """Test that remesh_keeping_periodicity_for_fem raises InputMeshNotPeriodicError."""
    if USE_MMG:
        edge_length_gradient = 1.05
        with pytest.raises(
            InputMeshNotPeriodicError,
            match="Input mesh is not periodic",
        ):
            remesh_keeping_periodicity_for_fem(
                non_periodic_mesh,
                hgrad=edge_length_gradient,
            )
