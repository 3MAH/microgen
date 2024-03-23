"""Tests for the remesh module."""

import subprocess
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

from .constants import (
    PERIODIC_BOX_FACES,
    PERIODIC_BOX_FACES_WITH_NB_NODES,
    PERIODIC_BOX_NODES,
)

_MESH_DIM = 3
_BOUNDARY_DIM = 2

USE_MMG = False
try:
    subprocess.run(["mmg3d_O3", "-h"], check=True)
    USE_MMG = True
except (subprocess.CalledProcessError, FileNotFoundError):
    warnings.warn("MMG will not be used in these tests")


@pytest.fixture(name="box_mesh", scope="function")
def fixture_box_mesh() -> BoxMesh:
    """Create a box mesh."""

    elements_dict = {pv.CellType.TETRA: PERIODIC_BOX_FACES}

    return BoxMesh(PERIODIC_BOX_NODES, elements_dict)


@pytest.fixture(name="gyroid_mesh", scope="function")
def fixture_gyroid_mesh() -> pv.UnstructuredGrid:
    """Create a gyroid mesh."""
    return pv.UnstructuredGrid(
        Tpms(surface_function=gyroid, offset=1.0).generateVtk(type_part="sheet")
    )


@pytest.fixture(name="non_periodic_mesh", scope="function")
def fixture_non_periodic_mesh() -> pv.UnstructuredGrid:
    """Create a non-periodic mesh."""
    elements = PERIODIC_BOX_FACES_WITH_NB_NODES
    cell_types = np.full(elements.shape[0], pv.CellType.TETRA, dtype=np.uint8)
    return pv.UnstructuredGrid(elements, cell_types, PERIODIC_BOX_NODES)


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
    """Test that remeshing a periodic mesh maintains periodicity."""
    # Arrange
    input_mesh = request.getfixturevalue(shape)
    # Act
    if USE_MMG:
        edge_length_gradient = 1.05
        remeshed_shape = remesh_keeping_periodicity_for_fem(
            input_mesh, hgrad=edge_length_gradient
        )

        if isinstance(input_mesh, BoxMesh):
            input_mesh = input_mesh.to_pyvista()
            remeshed_shape = remeshed_shape.to_pyvista()
        # Assert
        assert is_periodic(input_mesh.points) and is_periodic(remeshed_shape.points)


def test_given_non_periodic_mesh_remesh_must_raise_inputmeshnotperiodicerror(
    non_periodic_mesh: pv.UnstructuredGrid,
) -> None:
    """Test that remeshing a non-periodic mesh raises an error."""
    if USE_MMG:
        with pytest.raises(
            InputMeshNotPeriodicError, match="Input mesh is not periodic"
        ):
            edge_length_gradient = 1.05
            remesh_keeping_periodicity_for_fem(
                non_periodic_mesh, hgrad=edge_length_gradient
            )
