import pytest
import numpy as np
import numpy.typing as npt
from microgen import BoxMesh, Rve
import pyvista as pv


@pytest.fixture(scope='session')
def box_mesh_points() -> npt.NDArray[np.float_]:
    points = np.array([
        [1., 1., 1.],
        [1., 0.5, 1.],
        [1., 0., 1.],
        [1., 1., 0.5],
        [1., 0.5, 0.5],
        [1., 0., 0.5],
        [1., 1., 0.],
        [1., 0.5, 0.],
        [1., 0., 0.],
        [0.5, 1., 1.],
        [0.5, 0.5, 1.],
        [0.5, 0., 1.],
        [0.5, 1., 0.5],
        [0.5, 0.5, 0.5],
        [0.5, 0., 0.5],
        [0.5, 1., 0.],
        [0.5, 0.5, 0.],
        [0.5, 0., 0.],
        [0., 1., 1.],
        [0., 0.5, 1.],
        [0., 0., 1.],
        [0., 1., 0.5],
        [0., 0.5, 0.5],
        [0., 0., 0.5],
        [0., 1., 0.],
        [0., 0.5, 0.],
        [0., 0., 0.]
    ])
    return points


@pytest.fixture(scope="session")
def box_mesh(box_mesh_points) -> pv.UnstructuredGrid:
    points = box_mesh_points
    point_cloud = pv.PolyData(points)
    grid = point_cloud.delaunay_3d(offset=100.)

    return grid


@pytest.fixture(scope='session')
def box_boxMesh_elements() -> dict[int, npt.NDArray[np.float_]]:
    elements_dict = {pv.CellType.TETRA: np.array([[3, 4, 1, 9],
                                                  [3, 1, 0, 9],
                                                  [14, 16, 8, 17],
                                                  [6, 4, 3, 12],
                                                  [1, 9, 4, 10],
                                                  [4, 5, 2, 10],
                                                  [4, 2, 1, 10],
                                                  [4, 9, 3, 12],
                                                  [12, 10, 9, 18],
                                                  [10, 12, 4, 13],
                                                  [6, 7, 4, 12],
                                                  [5, 11, 10, 13],
                                                  [2, 10, 5, 11],
                                                  [8, 5, 7, 13],
                                                  [10, 9, 4, 12],
                                                  [4, 12, 7, 13],
                                                  [7, 5, 4, 13],
                                                  [5, 10, 4, 13],
                                                  [13, 12, 7, 15],
                                                  [12, 6, 7, 15],
                                                  [11, 13, 5, 14],
                                                  [5, 13, 8, 14],
                                                  [14, 13, 8, 16],
                                                  [13, 7, 8, 16],
                                                  [15, 7, 13, 16],
                                                  [12, 13, 10, 18],
                                                  [15, 13, 12, 21],
                                                  [10, 18, 13, 19],
                                                  [14, 11, 13, 19],
                                                  [13, 11, 10, 19],
                                                  [13, 18, 12, 21],
                                                  [14, 20, 19, 22],
                                                  [19, 21, 13, 22],
                                                  [19, 14, 11, 20],
                                                  [15, 16, 13, 21],
                                                  [21, 15, 16, 24],
                                                  [17, 14, 16, 22],
                                                  [19, 18, 13, 21],
                                                  [13, 21, 16, 22],
                                                  [16, 14, 13, 22],
                                                  [14, 19, 13, 22],
                                                  [20, 22, 14, 23],
                                                  [14, 22, 17, 23],
                                                  [22, 21, 16, 24],
                                                  [17, 22, 16, 25],
                                                  [22, 24, 16, 25],
                                                  [23, 22, 17, 25],
                                                  [23, 25, 17, 26]])}

    return elements_dict


@pytest.fixture(scope='session')
def box_BoxMesh(box_mesh_points, box_boxMesh_elements) -> BoxMesh:
    mesh = BoxMesh(box_mesh_points, box_boxMesh_elements)

    return mesh


def test_given_box_mesh_construct_must_find_center_corners_edges_faces_node_sets(box_BoxMesh) -> None:
    rve = Rve(dim_x=1.0, dim_y=1.0, dim_z=1.0, center=(0.5, 0.5, 0.5))
    target_center: npt.NDArray[np.float_] = np.array([0.5, 0.5, 0.5])
    target_corners: list[npt.NDArray[int]] = [np.array([26]),
                                              np.array([24]),
                                              np.array([8]),
                                              np.array([6]),
                                              np.array([20]),
                                              np.array([18]),
                                              np.array([2]),
                                              np.array([0])]
    target_edges: list[npt.NDArray[int]] = [np.array([23]),
                                            np.array([5]),
                                            np.array([3]),
                                            np.array([21]),
                                            np.array([25]),
                                            np.array([7]),
                                            np.array([1]),
                                            np.array([19]),
                                            np.array([17]),
                                            np.array([15]),
                                            np.array([9]),
                                            np.array([11])]
    target_faces: list[npt.NDArray[int]] = [np.array([22]), np.array([4]), np.array([14]), np.array([12]), np.array([16]), np.array([10])]

    box_BoxMesh.construct(rve)

    assert (box_BoxMesh.center == target_center).all() and box_BoxMesh.corners == target_corners and box_BoxMesh.edges == target_edges and box_BoxMesh.faces == target_faces
