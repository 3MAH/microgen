"""Tests for the is_periodic function."""

import numpy as np
import numpy.typing as npt
import pytest
import pyvista as pv

from microgen import is_periodic

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/
# ruff: noqa: E501 line-too-long https://docs.astral.sh/ruff/rules/line-too-long/


def _periodic_box_nodes() -> npt.NDArray[np.float64]:
    return np.array(
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


def _one_shifted_node_box_nodes() -> npt.NDArray[np.float64]:
    return np.array(
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


def _box_elements_same_number_of_nodes() -> npt.NDArray[np.int64]:
    return np.array(
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


def _one_extra_node_box_nodes() -> npt.NDArray[np.float64]:
    return np.array(
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
            [0.3, 0.2, -0.5],
        ],
    )


@pytest.fixture(name="periodic_box")
def fixture_periodic_box() -> pv.UnstructuredGrid:
    """Return a periodic box mesh."""
    nodes = _periodic_box_nodes()
    elements = _box_elements_same_number_of_nodes()
    cell_types = np.full(elements.shape[0], pv.CellType.TETRA, dtype=np.uint8)
    return pv.UnstructuredGrid(elements, cell_types, nodes)


@pytest.fixture(name="non_periodic_box_1_extra_node")
def fixture_non_periodic_box_1_extra_node() -> pv.UnstructuredGrid:
    """Return a non-periodic box mesh with an extra node."""
    points = _one_extra_node_box_nodes()
    point_cloud = pv.PolyData(points)
    return point_cloud.delaunay_3d(offset=100.0)


@pytest.fixture(name="non_periodic_box_shifted_node")
def fixture_non_periodic_box_shifted_node() -> pv.UnstructuredGrid:
    """Return a non-periodic box mesh with a shifted node but no extra node."""
    nodes = _one_shifted_node_box_nodes()
    elements = _box_elements_same_number_of_nodes()
    cell_types = np.full(elements.shape[0], pv.CellType.TETRA, dtype=np.uint8)
    return pv.UnstructuredGrid(elements, cell_types, nodes)


def test_given_periodic_box_is_periodic_must_return_true(
    periodic_box: pv.UnstructuredGrid,
) -> None:
    """is_periodic must return True for a periodic box."""
    crd = periodic_box.points

    assert is_periodic(crd)


def test_given_non_periodic_box_with_an_extra_node_is_periodic_must_return_false(
    non_periodic_box_1_extra_node: pv.UnstructuredGrid,
) -> None:
    """is_periodic must return False for a non-periodic box with an extra node."""
    crd = non_periodic_box_1_extra_node.points

    assert not is_periodic(crd)


def test_given_non_periodic_box_with_a_shifted_node_but_no_extra_node_is_periodic_must_return_false(
    non_periodic_box_shifted_node: pv.UnstructuredGrid,
) -> None:
    """is_periodic must return False for a non-periodic box with a shifted node."""
    crd = non_periodic_box_shifted_node.points

    assert not is_periodic(crd)
