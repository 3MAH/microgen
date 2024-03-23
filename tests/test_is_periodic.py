"""Tests for the is_periodic function."""

import numpy as np
import numpy.typing as npt
import pytest
import pyvista as pv

from microgen import is_periodic

from .constants import PERIODIC_BOX_FACES_WITH_NB_NODES, PERIODIC_BOX_NODES


def _one_shifted_node_box_nodes() -> npt.NDArray[np.float_]:
    box_nodes = PERIODIC_BOX_NODES.copy()
    box_nodes[-4, 1] = 0.1
    return box_nodes


def _one_extra_node_box_nodes() -> npt.NDArray[np.float_]:
    return np.vstack((PERIODIC_BOX_NODES, [0.3, 0.2, -0.5]))


@pytest.fixture(name="periodic_box", scope="function")
def fixture_periodic_box() -> pv.UnstructuredGrid:
    """Create a periodic box mesh."""
    cell_types = np.full(
        PERIODIC_BOX_FACES_WITH_NB_NODES.shape[0],
        pv.CellType.TETRA,
        dtype=np.uint8,
    )
    return pv.UnstructuredGrid(
        PERIODIC_BOX_FACES_WITH_NB_NODES,
        cell_types,
        PERIODIC_BOX_NODES,
    )


@pytest.fixture(name="non_periodic_box_1_extra_node", scope="function")
def fixture_non_periodic_box_1_extra_node() -> pv.UnstructuredGrid:
    """Create a non-periodic box mesh with an extra node."""
    points = _one_extra_node_box_nodes()
    point_cloud = pv.PolyData(points)
    return point_cloud.delaunay_3d(offset=100.0)


@pytest.fixture(name="non_periodic_box_shifted_node", scope="function")
def fixture_non_periodic_box_shifted_node() -> pv.UnstructuredGrid:
    """Create a non-periodic box mesh with a shifted node."""
    nodes = _one_shifted_node_box_nodes()
    cell_types = np.full(
        PERIODIC_BOX_FACES_WITH_NB_NODES.shape[0],
        pv.CellType.TETRA,
        dtype=np.uint8,
    )
    return pv.UnstructuredGrid(PERIODIC_BOX_FACES_WITH_NB_NODES, cell_types, nodes)


def test_given_periodic_box_is_periodic_must_return_true(
    periodic_box: pv.UnstructuredGrid,
) -> None:
    """Test if a periodic box is correctly identified as periodic."""
    crd = periodic_box.points

    assert is_periodic(crd)


def test_given_non_periodic_box_with_an_extra_node_is_periodic_must_return_false(
    non_periodic_box_1_extra_node: pv.UnstructuredGrid,
) -> None:
    """Test if a non-periodic box with an extra node is correctly identified as non-periodic."""
    crd = non_periodic_box_1_extra_node.points

    assert not is_periodic(crd)


def test_given_non_periodic_box_with_a_shifted_node_but_no_extra_node_is_periodic_must_return_false(
    non_periodic_box_shifted_node: pv.UnstructuredGrid,
) -> None:
    """Test if a non-periodic box with a shifted node but no extra node
    is correctly identified as non-periodic."""
    crd = non_periodic_box_shifted_node.points

    assert not is_periodic(crd)
