import pytest
import numpy as np
import numpy.typing as npt
from microgen import Box, PhaseMesh, meshPeriodic
import pyvista as pv
from pathlib import Path

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

def test_given_simple_periodic_pyvista_unstructured_grid_box_mesh_phaseMesh_from_pyvista_must_return_true(box_mesh_points, box_mesh) -> None:
    mesh = PhaseMesh.from_pyvista(box_mesh)

    assert mesh.nodes == box_mesh_points and mesh.elements == box_mesh.cells_dict

