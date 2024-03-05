import math as m

import numpy as np
import numpy.typing as npt
import pytest
import pyvista as pv
from microgen import BoxMesh, Rve


def _box_mesh_points() -> npt.NDArray[np.float_]:
    points = np.array(
        [
            [1.0, 1.0, 1.0],
            [1.0, 0.5, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 0.5],
            [1.0, 0.5, 0.5],
            [1.0, 0.0, 0.5],
            [1.0, 1.0, 0.0],
            [1.0, 0.5, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 1.0],
            [0.5, 0.5, 1.0],
            [0.5, 0.0, 1.0],
            [0.5, 1.0, 0.5],
            [0.5, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.0],
            [0.0, 1.0, 1.0],
            [0.0, 0.5, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.5],
            [0.0, 0.5, 0.5],
            [0.0, 0.0, 0.5],
            [0.0, 1.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )
    return points


def _box_mesh_elements() -> dict[pv.CellType, npt.NDArray[np.int_]]:
    elements_dict = {
        pv.CellType.TETRA: np.array(
            [
                [3, 4, 1, 9],
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
                [23, 25, 17, 26],
            ]
        )
    }

    return elements_dict


@pytest.fixture(name="default_box_mesh", scope="function")
def fixture_default_box_mesh() -> BoxMesh:
    mesh = BoxMesh(_box_mesh_points(), _box_mesh_elements())

    return mesh


def _default_rve() -> Rve:
    return Rve(dim_x=1.0, dim_y=1.0, dim_z=1.0, center=(0.5, 0.5, 0.5))


def _check_triangle_on_boundary(
    surface_mesh: pv.PolyData, triangle_index: int, rve: Rve
):
    triangle = surface_mesh.get_cell(triangle_index)
    triangle_nodes_coords = triangle.points.tolist()
    rve_boundaries = [
        (rve.x_min, rve.x_max),
        (rve.y_min, rve.y_max),
        (rve.z_min, rve.z_max),
    ]
    for i, rve_axis_min_max in enumerate(rve_boundaries):
        for rve_axis_boundary in rve_axis_min_max:
            is_triangle_on_boundary = all(
                m.isclose(node[i], rve_axis_boundary) for node in triangle_nodes_coords
            )
            if is_triangle_on_boundary:
                return True
    return False


def test_given_box_mesh__construct_must_find_center_corners_edges_faces_node_sets(
        default_box_mesh
) -> None:
    target_center: npt.NDArray[np.float_] = np.array([0.5, 0.5, 0.5])
    target_corners: list[int] = [0, 2, 6, 8, 18, 20, 24, 26]
    target_edges: list[int] = [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25]
    target_faces: list[int] = [4, 10, 12, 14, 16, 22]

    box_mesh_corners = sorted(np.concatenate(default_box_mesh.corners, axis=0).tolist())
    box_mesh_edges = sorted(np.concatenate(default_box_mesh.edges, axis=0).tolist())
    box_mesh_faces = sorted(np.concatenate(default_box_mesh.faces, axis=0).tolist())

    assert (
        (default_box_mesh.center == target_center).all()
        and box_mesh_corners == target_corners
        and box_mesh_edges == target_edges
        and box_mesh_faces == target_faces
    )


def test_given_box_mesh_rve_property_must_build_correct_rve(default_box_mesh) -> None:
    target_rve = _default_rve()
    test_rve = default_box_mesh.rve

    assert (
        test_rve.center == target_rve.center
        and test_rve.dim_x == target_rve.dim_x
        and test_rve.dim_y == target_rve.dim_y
        and test_rve.dim_z == target_rve.dim_z
    )


def test_given_box_box_mesh_boundary_elements_must_find_boundary_surface_elements(
        default_box_mesh
) -> None:
    expected_number_of_cells = 48
    rve = _default_rve()
    boundary, boundary_cells_index = default_box_mesh.boundary_elements(rve)
    bool_check_triangle_on_boundary_list = []
    for triangle_index in boundary_cells_index:
        bool_check_triangle_on_boundary_list.append(
            _check_triangle_on_boundary(boundary, triangle_index, rve)
        )

    assert boundary.n_faces == expected_number_of_cells and all(
        bool_check_triangle_on_boundary_list
    )
