"""Tests for the SingleMesh class."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pytest
import pyvista as pv

from microgen import SingleMesh, check_if_only_linear_tetrahedral
from microgen.single_mesh import NotOnlyLinearTetrahedraError

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/
# ruff: noqa: E501 line-too-long https://docs.astral.sh/ruff/rules/line-too-long/


def are_elements_equal(
    dict1: dict[pv.CellType, npt.NDArray[np.int_]],
    dict2: dict[pv.CellType, npt.NDArray[np.int_]],
) -> bool:
    """Return whether two dictionaries of arrays are equal."""
    if dict1.keys() != dict2.keys():
        return False
    return all(np.array_equal(dict1[key], dict2[key]) for key in dict1)


def _box_mesh_points() -> npt.NDArray[np.float64]:
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


@pytest.fixture(name="simple_grid")
def fixture_simple_grid() -> pv.UnstructuredGrid:
    """Return a simple grid."""
    points = _box_mesh_points()
    point_cloud = pv.PolyData(points)
    return point_cloud.delaunay_3d(offset=100.0)


def _box_single_mesh_nodes() -> npt.NDArray[np.float64]:
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


def _box_single_mesh_elements() -> dict[pv.CellType, npt.NDArray[np.int_]]:
    return {
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


@pytest.fixture(name="box_mesh")
def fixture_box_mesh() -> SingleMesh:
    """Return a box mesh."""
    return SingleMesh(_box_single_mesh_nodes(), _box_single_mesh_elements())


def _linear_1d_mesh() -> pv.UnstructuredGrid:
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    cells = [2, 0, 1]
    celltypes = [pv.CellType.LINE]
    return pv.UnstructuredGrid(cells, celltypes, points)


def _quadratic_1d_mesh() -> pv.UnstructuredGrid:
    points = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]])
    cells = [3, 0, 1, 2]
    celltypes = [pv.CellType.QUADRATIC_EDGE]
    return pv.UnstructuredGrid(cells, celltypes, points)


def _linear_2d_quad_mesh() -> pv.UnstructuredGrid:
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
    )
    cells = [4, 0, 1, 2, 3]
    celltypes = [pv.CellType.QUAD]
    return pv.UnstructuredGrid(cells, celltypes, points)


def _quadratic_2d_quad_mesh() -> pv.UnstructuredGrid:
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.0, 0.0],
            [1.0, 0.5, 0.0],
            [0.5, 1.0, 0.0],
            [0.0, 0.5, 0.0],
        ],
    )
    cells = [8, 0, 1, 2, 3, 4, 5, 6, 7]
    celltypes = [pv.CellType.QUADRATIC_QUAD]
    return pv.UnstructuredGrid(cells, celltypes, points)


def _linear_2d_triangle_mesh() -> pv.UnstructuredGrid:
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    cells = [3, 0, 1, 2]
    celltypes = [pv.CellType.TRIANGLE]
    return pv.UnstructuredGrid(cells, celltypes, points)


def _quadratic_2d_triangle_mesh() -> pv.UnstructuredGrid:
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.0],
        ],
    )
    cells = [6, 0, 1, 2, 3, 4, 5]
    celltypes = [pv.CellType.QUADRATIC_TRIANGLE]
    return pv.UnstructuredGrid(cells, celltypes, points)


def _linear_3d_hex_mesh() -> pv.UnstructuredGrid:
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
    )
    cells = [8, 0, 1, 2, 3, 4, 5, 6, 7]
    celltypes = [pv.CellType.HEXAHEDRON]
    return pv.UnstructuredGrid(cells, celltypes, points)


def _quadratic_3d_hex_mesh() -> pv.UnstructuredGrid:
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
            [0.5, 0.0, 0.0],
            [1.0, 0.5, 0.0],
            [0.5, 1.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.5, 0.0, 1.0],
            [1.0, 0.5, 1.0],
            [0.5, 1.0, 1.0],
            [0.0, 0.5, 1.0],
            [0.0, 0.0, 0.5],
            [1.0, 0.0, 0.5],
            [1.0, 1.0, 0.5],
            [0.0, 1.0, 0.5],
        ],
    )
    cells = [20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
    celltypes = [pv.CellType.QUADRATIC_HEXAHEDRON]
    return pv.UnstructuredGrid(cells, celltypes, points)


def _quadratic_3d_tet_mesh() -> pv.UnstructuredGrid:
    points = np.array(
        [
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, -1.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
            [-1.0, 0.0, 0.0],
        ],
    )
    cells = [10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    celltypes = [pv.CellType.QUADRATIC_TETRA]
    return pv.UnstructuredGrid(cells, celltypes, points)


def _linear_3d_wedge_mesh() -> pv.UnstructuredGrid:
    points = np.array(
        [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.5, 0.5],
        ],
    )
    cells = [6, 0, 1, 2, 3, 4, 5]
    celltypes = [pv.CellType.WEDGE]
    return pv.UnstructuredGrid(cells, celltypes, points)


def _linear_3d_pyramid_mesh() -> pv.UnstructuredGrid:
    points = np.array(
        [
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
            [0.0, 0.0, 1.60803807],
        ],
    )
    cells = [5, 0, 1, 2, 3, 4]
    celltypes = [pv.CellType.PYRAMID]
    return pv.UnstructuredGrid(cells, celltypes, points)


@pytest.fixture(name="sample_1d_mesh_list")
def fixture_sample_1d_mesh_list() -> list[pv.UnstructuredGrid]:
    """Return a list of 1D meshes with linear and quadratic elements."""
    return [_linear_1d_mesh(), _quadratic_1d_mesh()]


@pytest.fixture(name="sample_2d_mesh_list")
def fixture_sample_2d_mesh_list() -> list[pv.UnstructuredGrid]:
    """Return a list of 2D meshes with linear and quadratic elements."""
    return [
        _linear_2d_quad_mesh(),
        _quadratic_2d_quad_mesh(),
        _linear_2d_triangle_mesh(),
        _quadratic_2d_triangle_mesh(),
    ]


@pytest.fixture(name="sample_3d_non_linear_tet_mesh_list")
def fixture_sample_3d_non_linear_tet_mesh_list() -> list[pv.UnstructuredGrid]:
    """Return a list of 3D meshes with non-linear tetrahedral elements."""
    return [
        _quadratic_3d_tet_mesh(),
        _linear_3d_hex_mesh(),
        _quadratic_3d_hex_mesh(),
        _linear_3d_wedge_mesh(),
        _linear_3d_pyramid_mesh(),
    ]


def test_check_single_mesh_consistency_when_importing_mesh_from_pyvista(
    simple_grid: pv.UnstructuredGrid,
) -> None:
    """Test that the mesh is imported correctly."""
    mesh = SingleMesh.from_pyvista(simple_grid)

    assert mesh.nodes_coords.all() == simple_grid.points.all()
    assert are_elements_equal(mesh.elements, simple_grid.cells_dict)


def test_check_single_mesh_consistency_when_exporting_mesh_to_pyvista(
    box_mesh: SingleMesh,
) -> None:
    """Test that the mesh is exported correctly."""
    grid = box_mesh.to_pyvista()

    assert are_elements_equal(grid.cells_dict, box_mesh.elements)


def test_single_mesh_consistency_when_extracting_surface_mesh(
    simple_grid: pv.UnstructuredGrid,
) -> None:
    """Test that the surface mesh is extracted correctly."""
    mesh = SingleMesh.from_pyvista(simple_grid)

    target_n_cells = 24
    target_face_connectivity_array = np.array(
        [
            [3, 0, 1, 2],
            [3, 0, 3, 1],
            [3, 0, 4, 5],
            [3, 0, 2, 4],
            [3, 0, 5, 6],
            [3, 0, 6, 3],
            [3, 7, 8, 9],
            [3, 7, 9, 10],
            [3, 7, 11, 8],
            [3, 7, 12, 11],
            [3, 7, 10, 13],
            [3, 7, 13, 12],
            [3, 9, 8, 1],
            [3, 9, 1, 3],
            [3, 9, 6, 10],
            [3, 9, 3, 6],
            [3, 4, 2, 11],
            [3, 4, 11, 12],
            [3, 4, 13, 5],
            [3, 4, 12, 13],
            [3, 1, 8, 11],
            [3, 1, 11, 2],
            [3, 6, 5, 13],
            [3, 6, 13, 10],
        ],
    ).reshape(-1)

    surf = mesh.surface

    assert surf.n_cells == target_n_cells
    assert surf.faces.all() == target_face_connectivity_array.all()


def test_given_sample_1d_mesh__check_if_only_linear_tetrahedral_must_raise_1d_warning(
    sample_1d_mesh_list: list[pv.UnstructuredGrid],
) -> None:
    """Test that the function raises a warning when the mesh contains 1D elements."""
    warning_message = (
        "1D elements are present in the PyVista UnstructuredGrid. They will be ignored."
    )
    with pytest.warns(UserWarning, match=warning_message):
        for mesh in sample_1d_mesh_list:
            check_if_only_linear_tetrahedral(mesh)


def test_given_sample_2d_mesh__check_if_only_linear_tetrahedral_must_raise_2d_warning(
    sample_2d_mesh_list: list[pv.UnstructuredGrid],
) -> None:
    """Test that the function raises a warning when the mesh contains 2D elements."""
    warning_message = (
        "2D elements are present in the PyVista UnstructuredGrid. They will be ignored."
    )
    with pytest.warns(UserWarning, match=warning_message):
        for mesh in sample_2d_mesh_list:
            check_if_only_linear_tetrahedral(mesh)


def test_given_sample_3d_mesh__check_if_only_linear_tetrahedral_must_raise_found_non_linear_tet_elements_error(
    sample_3d_non_linear_tet_mesh_list: list[pv.UnstructuredGrid],
) -> None:
    """Test that the function raises an error when the mesh contains non-linear tetrahedra."""
    error_message_snippet = "Mesh contains elements other than linear tetrahedra."
    for mesh in sample_3d_non_linear_tet_mesh_list:
        with pytest.raises(NotOnlyLinearTetrahedraError, match=error_message_snippet):
            check_if_only_linear_tetrahedral(mesh)


def test_given_linear_tets_only_mesh__check_if_only_linear_tetrahedral_must_not_raise_any_error(
    simple_grid: pv.UnstructuredGrid,
) -> None:
    """Test that the function does not raise an error when the mesh contains only linear tetrahedra."""
    check_if_only_linear_tetrahedral(simple_grid)
