"""Tests for the BoxMesh class."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pytest
import pyvista as pv

from microgen import BoxMesh, Rve

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/
# ruff: noqa: E501 line-too-long https://docs.astral.sh/ruff/rules/line-too-long/


def _box_mesh_points() -> npt.NDArray[np.float64]:
    return np.array(
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
        ],
    )


def _box_mesh_elements() -> dict[pv.CellType, npt.NDArray[np.int_]]:
    return {
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
            ],
        ),
    }


@pytest.fixture()
def default_box_mesh() -> tuple[BoxMesh, Rve]:
    """Return a default box mesh and its corresponding RVE."""
    mesh = BoxMesh(_box_mesh_points(), _box_mesh_elements())
    rve = Rve(dim=(1, 1, 1), center=(0.5, 0.5, 0.5))
    return (mesh, rve)


@pytest.fixture()
def expanded_box_mesh_x() -> tuple[BoxMesh, Rve]:
    """Return a box mesh expanded in the x direction and its corresponding RVE."""
    expanded_mesh = _box_mesh_points()
    expanded_mesh[:, 0] *= 3.0
    mesh = BoxMesh(expanded_mesh, _box_mesh_elements())
    rve = Rve(dim=(3, 1, 1), center=(1.5, 0.5, 0.5))
    return (mesh, rve)


@pytest.fixture()
def expanded_box_mesh_y() -> tuple[BoxMesh, Rve]:
    """Return a box mesh expanded in the y direction and its corresponding RVE."""
    expanded_mesh = _box_mesh_points()
    expanded_mesh[:, 1] *= 3.0
    mesh = BoxMesh(expanded_mesh, _box_mesh_elements())
    rve = Rve(dim=(1, 3, 1), center=(0.5, 1.5, 0.5))
    return (mesh, rve)


@pytest.fixture()
def expanded_box_mesh_z() -> tuple[BoxMesh, Rve]:
    """Return a box mesh expanded in the z direction and its corresponding RVE."""
    expanded_mesh = _box_mesh_points()
    expanded_mesh[:, 2] *= 3.0
    mesh = BoxMesh(expanded_mesh, _box_mesh_elements())
    rve = Rve(dim=(1, 1, 3), center=(0.5, 0.5, 1.5))
    return (mesh, rve)


def _check_triangle_on_boundary(
    surface_mesh: pv.PolyData,
    triangle_index: int,
    rve: Rve,
) -> bool:
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
                np.isclose(node[i], rve_axis_boundary) for node in triangle_nodes_coords
            )
            if is_triangle_on_boundary:
                return True
    return False


@pytest.mark.parametrize(
    "box_mesh",
    [
        "default_box_mesh",
        "expanded_box_mesh_x",
        "expanded_box_mesh_y",
        "expanded_box_mesh_z",
    ],
)
def test_given_box_mesh_construct_must_find_center_corners_edges_faces_node_sets(
    box_mesh: BoxMesh,
    request: pytest.FixtureRequest,
) -> None:
    """Test that the BoxMesh object is correctly constructed."""
    box_mesh_to_test, rve = request.getfixturevalue(box_mesh)
    target_center: npt.NDArray[np.float64] = rve.center
    target_corners: list[int] = [0, 2, 6, 8, 18, 20, 24, 26]
    target_edges: list[int] = [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25]
    target_faces: list[int] = [4, 10, 12, 14, 16, 22]

    box_mesh_corners = sorted(np.concatenate(list(box_mesh_to_test.corners.values())))
    box_mesh_edges = sorted(np.concatenate(list(box_mesh_to_test.edges.values())))
    box_mesh_faces = sorted(np.concatenate(list(box_mesh_to_test.faces.values())))

    assert (box_mesh_to_test.center == target_center).all()
    assert box_mesh_corners == target_corners
    assert box_mesh_edges == target_edges
    assert box_mesh_faces == target_faces


@pytest.mark.parametrize(
    "box_mesh",
    [
        "default_box_mesh",
        "expanded_box_mesh_x",
        "expanded_box_mesh_y",
        "expanded_box_mesh_z",
    ],
)
def test_given_box_mesh_rve_property_must_build_correct_rve(
    box_mesh: BoxMesh,
    request: pytest.FixtureRequest,
) -> None:
    """Test that the BoxMesh object builds the correct RVE."""
    box_mesh_to_test, target_rve = request.getfixturevalue(box_mesh)
    test_rve = box_mesh_to_test.rve

    assert np.all(test_rve.center == target_rve.center)
    assert np.all(test_rve.dim == target_rve.dim)


@pytest.mark.parametrize(
    "box_mesh",
    [
        "default_box_mesh",
        "expanded_box_mesh_x",
        "expanded_box_mesh_y",
        "expanded_box_mesh_z",
    ],
)
def test_given_box_box_mesh_boundary_elements_must_find_boundary_surface_elements(
    box_mesh: BoxMesh,
    request: pytest.FixtureRequest,
) -> None:
    """Test that the BoxMesh object finds the correct boundary elements."""
    box_mesh_to_test, rve = request.getfixturevalue(box_mesh)
    expected_number_of_cells = 48
    boundary, boundary_cells_index = box_mesh_to_test.boundary_elements(rve)
    bool_check_triangle_on_boundary_list = [
        _check_triangle_on_boundary(boundary, triangle_index, rve)
        for triangle_index in boundary_cells_index
    ]
    assert boundary.n_cells == expected_number_of_cells
    assert all(bool_check_triangle_on_boundary_list)


@pytest.mark.parametrize("k_neighbours", [1, 2, 3, 4])
def test_closest_points_on_boundaries_periodic_mesh_must_have_only_one_neighbour(
    k_neighbours: int,
) -> None:
    """Test that the closest points on the boundaries of a periodic mesh have only one neighbour."""
    box_mesh_to_test = BoxMesh(_box_mesh_points(), _box_mesh_elements())
    closest_pts = box_mesh_to_test.closest_points_on_boundaries(
        k_neighbours=k_neighbours,
    )

    number_of_closest_neighbours_to_test: list[int] = []
    number_closest_neighbours_truth_1_on_each_face_and_edge = [1] * 12
    for neighbours_for_all_nodes in closest_pts.values():
        neighbours_for_all_nodes_index = neighbours_for_all_nodes[0]

        number_of_closest_neighbours_to_test.extend(
            len(neighbours_for_each_nodes)
            for neighbours_for_each_nodes in neighbours_for_all_nodes_index
        )

    assert (
        number_of_closest_neighbours_to_test
        == number_closest_neighbours_truth_1_on_each_face_and_edge
    )


@pytest.mark.parametrize("k_neighbours", [1, 2, 3, 4])
def test_closest_points_on_perturbated_mesh_boundaries_must_have_good_number_of_neighbours(
    k_neighbours: int,
) -> None:
    """Test that the closest points on the boundaries of a perturbated mesh have the correct number of neighbours."""
    box_mesh_to_test = BoxMesh(_box_mesh_points(), _box_mesh_elements())
    perturbation_factor = 0.01

    xyz_min_values_of_rve = np.min(box_mesh_to_test.rve.min_point)
    xyz_max_values_of_rve = np.max(box_mesh_to_test.rve.max_point)

    rng = np.random.default_rng()

    # perturb faces and edges with non-zero perturbation so that multiple neighbours must be found
    for node_index_in_faces in box_mesh_to_test.faces.values():
        box_mesh_to_test.nodes_coords[node_index_in_faces[0]] += perturbation_factor * (
            1 - rng.uniform(0, 1)
        )

    for node_index_in_edges in box_mesh_to_test.edges.values():
        box_mesh_to_test.nodes_coords[node_index_in_edges[0]] += perturbation_factor * (
            1 - rng.uniform(0, 1)
        )

    box_mesh_to_test.nodes_coords = np.clip(
        box_mesh_to_test.nodes_coords,
        xyz_min_values_of_rve,
        xyz_max_values_of_rve,
    )

    closest_pts = box_mesh_to_test.closest_points_on_boundaries(
        k_neighbours=k_neighbours,
    )

    number_closest_neighbours_truth_k_neighbours_on_each_face_2_on_each_edge = [
        k_neighbours,
        k_neighbours,
        k_neighbours,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
    ]
    number_of_closest_neighbours_to_test: list[int] = []
    for neighbours_for_all_nodes in closest_pts.values():
        neighbours_for_all_nodes_index = neighbours_for_all_nodes[0]

        number_of_closest_neighbours_to_test.extend(
            len(neighbours_for_each_nodes)
            for neighbours_for_each_nodes in neighbours_for_all_nodes_index
        )

    assert (
        number_of_closest_neighbours_to_test
        == number_closest_neighbours_truth_k_neighbours_on_each_face_2_on_each_edge
    )
