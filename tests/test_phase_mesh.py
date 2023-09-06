import pytest
import numpy as np
import numpy.typing as npt
from microgen import PhaseMesh, phaseMesh
import pyvista as pv
import warnings


def compare_dict_with_arrays_as_values(dict1: dict[int, npt.NDArray[int]], dict2: dict[int, npt.NDArray[int]]) -> bool:
    """Return whether two dictionaries of arrays are equal"""
    if dict1.keys() != dict2.keys():
        return False
    return all(np.array_equal(dict1[key], dict2[key]) for key in dict1)


@pytest.fixture(scope="session")
def box_mesh_points() -> npt.NDArray[np.float_]:
    points = np.array([
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
        [0.0, 0.0, -0.5]
    ])
    return points


@pytest.fixture(scope="session")
def box_mesh(box_mesh_points) -> pv.UnstructuredGrid:
    points = box_mesh_points
    point_cloud = pv.PolyData(points)
    grid = point_cloud.delaunay_3d(offset=100.)

    return grid


@pytest.fixture(scope='session')
def box_phaseMesh_nodes() -> npt.NDArray[np.float_]:
    nodes_array = np.array([[0., 0., 0.],
                            [0.5, 0.5, 0.5],
                            [-0.5, -0.5, -0.5],
                            [-0.5, -0.5, 0.5],
                            [0.5, 0.5, -0.5],
                            [0.5, -0.5, 0.5],
                            [0.5, -0.5, -0.5],
                            [-0.5, 0.5, 0.5],
                            [-0.5, 0.5, -0.5],
                            [0., 0.5, 0.],
                            [0., -0.5, 0.],
                            [0.5, 0., 0.],
                            [-0.5, 0., 0.],
                            [0., 0., 0.5],
                            [0., 0., -0.5]])

    return nodes_array


@pytest.fixture(scope='session')
def box_phaseMesh_elements() -> dict[int, npt.NDArray[np.float_]]:
    elements_dict = {10: np.array([[11, 6, 0, 14],
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
                                   [9, 8, 4, 14]])}

    return elements_dict


@pytest.fixture(scope='session')
def box_phaseMesh(box_phaseMesh_nodes, box_phaseMesh_elements) -> PhaseMesh:
    box_phasemesh = PhaseMesh(box_phaseMesh_nodes, box_phaseMesh_elements)

    return box_phasemesh


@pytest.fixture(scope='session')
def linear_1d_mesh() -> pv.UnstructuredGrid:
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    cells = [2, 0, 1]
    celltypes = [pv.CellType.LINE]
    grid = pv.UnstructuredGrid(cells, celltypes, points)

    return grid


@pytest.fixture(scope='session')
def quadratic_1d_mesh() -> pv.UnstructuredGrid:
    points = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]])
    cells = [3, 0, 1, 2]
    celltypes = [pv.CellType.QUADRATIC_EDGE]
    grid = pv.UnstructuredGrid(cells, celltypes, points)

    return grid


@pytest.fixture(scope='session')
def linear_2d_quad_mesh() -> pv.UnstructuredGrid:
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    cells = [4, 0, 1, 2, 3]
    celltypes = [pv.CellType.QUAD]
    grid = pv.UnstructuredGrid(cells, celltypes, points)

    return grid


@pytest.fixture(scope='session')
def quadratic_2d_quad_mesh() -> pv.UnstructuredGrid:
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.5, 0.0],
         [0.5, 1.0, 0.0], [0.0, 0.5, 0.0]])
    cells = [8, 0, 1, 2, 3, 4, 5, 6, 7]
    celltypes = [pv.CellType.QUADRATIC_QUAD]
    grid = pv.UnstructuredGrid(cells, celltypes, points)

    return grid


@pytest.fixture(scope='session')
def linear_2d_triangle_mesh() -> pv.UnstructuredGrid:
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    cells = [3, 0, 1, 2]
    celltypes = [pv.CellType.TRIANGLE]
    grid = pv.UnstructuredGrid(cells, celltypes, points)

    return grid


@pytest.fixture(scope='session')
def quadratic_2d_triangle_mesh() -> pv.UnstructuredGrid:
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.0]])
    cells = [6, 0, 1, 2, 3, 4, 5]
    celltypes = [pv.CellType.QUADRATIC_TRIANGLE]
    grid = pv.UnstructuredGrid(cells, celltypes, points)

    return grid


@pytest.fixture(scope='session')
def linear_3d_hex_mesh() -> pv.UnstructuredGrid:
    points = np.array([[0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0],
                       [1.0, 1.0, 0.0],
                       [0.0, 1.0, 0.0],
                       [0.0, 0.0, 1.0],
                       [1.0, 0.0, 1.0],
                       [1.0, 1.0, 1.0],
                       [0.0, 1.0, 1.0]])
    cells = [8, 0, 1, 2, 3, 4, 5, 6, 7]
    celltypes = [pv.CellType.HEXAHEDRON]
    grid = pv.UnstructuredGrid(cells, celltypes, points)

    return grid


@pytest.fixture(scope='session')
def quadratic_3d_hex_mesh() -> pv.UnstructuredGrid:
    points = np.array([[0.0, 0.0, 0.0],
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
                       ])
    cells = [20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
    celltypes = [pv.CellType.QUADRATIC_HEXAHEDRON]
    grid = pv.UnstructuredGrid(cells, celltypes, points)

    return grid


@pytest.fixture(scope='session')
def quadratic_3d_tet_mesh() -> pv.UnstructuredGrid:
    points = np.array([[1.0, 1.0, 1.0],
                       [1.0, -1.0, -1.0],
                       [-1.0, 1.0, -1.0],
                       [-1.0, -1.0, 1.0],
                       [1.0, 0.0, 0.0],
                       [0.0, 0.0, -1.0],
                       [0.0, 1.0, 0.0],
                       [0.0, 0.0, 1.0],
                       [0.0, -1.0, 0.0],
                       [-1.0, 0.0, 0.0]
                       ])
    cells = [10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    celltypes = [pv.CellType.QUADRATIC_TETRA]
    grid = pv.UnstructuredGrid(cells, celltypes, points)

    return grid


@pytest.fixture(scope='session')
def linear_3d_wedge_mesh() -> pv.UnstructuredGrid:
    points = np.array([[0.0, 1.0, 0.0],
                       [0.0, 0.0, 0.0],
                       [0.0, 0.5, 0.5],
                       [1.0, 1.0, 0.0],
                       [1.0, 0.0, 0.0],
                       [1.0, 0.5, 0.5]])
    cells = [6, 0, 1, 2, 3, 4, 5]
    celltypes = [pv.CellType.WEDGE]
    grid = pv.UnstructuredGrid(cells, celltypes, points)

    return grid


@pytest.fixture(scope='session')
def linear_3d_pyramid_mesh() -> pv.UnstructuredGrid:
    points = np.array([[1.0, 1.0, 0.0],
                       [-1.0, 1.0, 0.0],
                       [-1.0, -1.0, 0.0],
                       [1.0, -1.0, 0.0],
                       [0.0, 0.0, 1.60803807]])
    cells = [5, 0, 1, 2, 3, 4]
    celltypes = [pv.CellType.PYRAMID]
    grid = pv.UnstructuredGrid(cells, celltypes, points)

    return grid


@pytest.fixture(scope='session')
def sample_1d_mesh_list(linear_1d_mesh: pv.UnstructuredGrid, quadratic_1d_mesh: pv.UnstructuredGrid) -> list[
    pv.UnstructuredGrid]:
    return [linear_1d_mesh, quadratic_1d_mesh]


@pytest.fixture(scope='session')
def sample_2d_mesh_list(linear_2d_quad_mesh: pv.UnstructuredGrid, quadratic_2d_quad_mesh, linear_2d_triangle_mesh,
                        quadratic_2d_triangle_mesh) -> list[pv.UnstructuredGrid]:
    return [linear_2d_quad_mesh, quadratic_2d_quad_mesh, linear_2d_triangle_mesh, quadratic_2d_triangle_mesh]


@pytest.fixture(scope='session')
def sample_3d_non_linear_tet_mesh_list(quadratic_3d_tet_mesh: pv.UnstructuredGrid,
                                       linear_3d_hex_mesh: pv.UnstructuredGrid,
                                       quadratic_3d_hex_mesh: pv.UnstructuredGrid,
                                       linear_3d_wedge_mesh: pv.UnstructuredGrid,
                                       linear_3d_pyramid_mesh: pv.UnstructuredGrid) -> list[pv.UnstructuredGrid]:
    return [quadratic_3d_tet_mesh, linear_3d_hex_mesh, quadratic_3d_hex_mesh, linear_3d_wedge_mesh, linear_3d_pyramid_mesh]


def test_given_simple_periodic_pyvista_unstructured_grid_box_mesh_phaseMesh_from_pyvista_must_return_the_same_mesh(
        box_mesh_points, box_mesh) -> None:
    mesh = PhaseMesh.from_pyvista(box_mesh)

    assert mesh.nodes.all() == box_mesh_points.all() and compare_dict_with_arrays_as_values(mesh.elements,
                                                                                            box_mesh.cells_dict)


def test_given_simple_periodic_box_phaseMesh_to_pyvista_must_return_the_same_mesh(box_phaseMesh) -> None:
    grid = box_phaseMesh.to_pyvista()

    assert compare_dict_with_arrays_as_values(grid.cells_dict, box_phaseMesh.elements)


def test_given_simple_periodic_pyvista_unstructured_grid_box_mesh_phaseMesh_surface_must_find_surface_triangles_connectivity_array_and_number(
        box_mesh) -> None:
    mesh = PhaseMesh.from_pyvista(box_mesh)

    target_n_cells = 24
    target_face_connectivity_array = np.array([3, 0, 1, 2, 3, 0, 3, 1, 3, 0, 4, 5, 3, 0, 2, 4, 3,
                                               0, 5, 6, 3, 0, 6, 3, 3, 7, 8, 9, 3, 7, 9, 10, 3, 7,
                                               11, 8, 3, 7, 12, 11, 3, 7, 10, 13, 3, 7, 13, 12, 3, 9, 8,
                                               1, 3, 9, 1, 3, 3, 9, 6, 10, 3, 9, 3, 6, 3, 4, 2, 11,
                                               3, 4, 11, 12, 3, 4, 13, 5, 3, 4, 12, 13, 3, 1, 8, 11, 3,
                                               1, 11, 2, 3, 6, 5, 13, 3, 6, 13, 10])

    surf = mesh.surface

    assert surf.n_faces == target_n_cells and surf.faces.all() == target_face_connectivity_array.all()

def test_given_sample_1d_mesh__check_if_only_linear_tetrahedral_must_raise_1d_warning(sample_1d_mesh_list) -> None:
    warning_message = "1D elements are present in the PyVista UnstructuredGrid. They will be ignored."
    with pytest.warns(UserWarning, match=warning_message):
        for mesh in sample_1d_mesh_list:
            phaseMesh._check_if_only_linear_tetrahedral(mesh)

def test_given_sample_2d_mesh__check_if_only_linear_tetrahedral_must_raise_2d_warning(sample_2d_mesh_list) -> None:
    warning_message = "2D elements are present in the PyVista UnstructuredGrid. They will be ignored. \nSurface elements shall be extracted automatically from the 3d mesh"
    with pytest.warns(UserWarning, match=warning_message):
        for mesh in sample_2d_mesh_list:
            phaseMesh._check_if_only_linear_tetrahedral(mesh)