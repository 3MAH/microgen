import pytest
from microgen import Box, meshPeriodic, Phase, Rve, remesh, Tpms, tpms, BasicGeometry
import cadquery as cq
import numpy as np
import numpy.typing as npt
import gmsh
from pathlib import Path

_MESH_DIM = 3
_BOUNDARY_DIM = 2

USE_MMG = False
try:
    subprocess.check_output("mmg3d_O3", stderr=subprocess.STDOUT)
    USE_MMG= True    
except subprocess.CalledProcessError:
    print(
        "mmg command did not work, check if it is installed or contact a developer"
    )
    USE_MMG= False 

def _get_mesh_nodes_coords(mesh_name: str) -> npt.NDArray[np.float_]:
    gmsh.initialize()
    gmsh.open(mesh_name)

    _, nodes_coords, _ = gmsh.model.mesh.getNodes()
    gmsh.finalize()

    nodes_coords = nodes_coords.reshape(-1, _MESH_DIM)

    return nodes_coords


def _is_periodic(nodes_coords, tol=1e-8, dim=_MESH_DIM) -> bool:
    # bounding box
    xmax = np.max(nodes_coords[:, 0])
    xmin = np.min(nodes_coords[:, 0])
    ymax = np.max(nodes_coords[:, 1])
    ymin = np.min(nodes_coords[:, 1])
    if dim == 3:
        zmax = np.max(nodes_coords[:, 2])
        zmin = np.min(nodes_coords[:, 2])

    # extract face nodes
    left = np.where(np.abs(nodes_coords[:, 0] - xmin) < tol)[0]
    right = np.where(np.abs(nodes_coords[:, 0] - xmax) < tol)[0]

    if dim > 1:
        bottom = np.where(np.abs(nodes_coords[:, 1] - ymin) < tol)[0]
        top = np.where(np.abs(nodes_coords[:, 1] - ymax) < tol)[0]

    if dim > 2:  # or dim == 3
        back = np.where(np.abs(nodes_coords[:, 2] - zmin) < tol)[0]
        front = np.where(np.abs(nodes_coords[:, 2] - zmax) < tol)[0]

        # sort adjacent faces to ensure node correspondance
    if nodes_coords.shape[1] == 2:  # 2D mesh
        left = left[np.argsort(nodes_coords[left, 1])]
        right = right[np.argsort(nodes_coords[right, 1])]
        if dim > 1:
            bottom = bottom[np.argsort(nodes_coords[bottom, 0])]
            top = top[np.argsort(nodes_coords[top, 0])]

    elif nodes_coords.shape[1] > 2:
        decimal_round = int(-np.log10(tol) - 1)
        left = left[
            np.lexsort(
                (nodes_coords[left, 1], nodes_coords[left, 2].round(decimal_round))
            )
        ]
        right = right[
            np.lexsort(
                (nodes_coords[right, 1], nodes_coords[right, 2].round(decimal_round))
            )
        ]
        if dim > 1:
            bottom = bottom[
                np.lexsort(
                    (
                        nodes_coords[bottom, 0],
                        nodes_coords[bottom, 2].round(decimal_round),
                    )
                )
            ]
            top = top[
                np.lexsort(
                    (nodes_coords[top, 0], nodes_coords[top, 2].round(decimal_round))
                )
            ]
        if dim > 2:
            back = back[
                np.lexsort(
                    (nodes_coords[back, 0], nodes_coords[back, 1].round(decimal_round))
                )
            ]
            front = front[
                np.lexsort(
                    (
                        nodes_coords[front, 0],
                        nodes_coords[front, 1].round(decimal_round),
                    )
                )
            ]

    # ==========================
    # test if mesh is periodic:
    # ==========================

    # test if same number of nodes in adjacent faces
    if len(left) != len(right):
        return False
    if dim > 1 and len(bottom) != len(top):
        return False
    if dim > 2 and (len(back) != len(front)):
        return False

    # check nodes position
    if (nodes_coords[right, 1:] - nodes_coords[left, 1:] > tol).any():
        return False
    if dim > 1 and (nodes_coords[top, ::2] - nodes_coords[bottom, ::2] > tol).any():
        return False
    if dim > 2 and (nodes_coords[front, :2] - nodes_coords[back, :2] > tol).any():
        return False

    return True


@pytest.fixture(scope="function")
def tmp_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    tmp_dir_name = "test_tmp_dir"
    return tmp_path_factory.mktemp(tmp_dir_name)


@pytest.fixture(scope="function")
def tmp_step_filename(tmp_dir: Path) -> str:
    return (tmp_dir / "shape.step").as_posix()


@pytest.fixture(scope="function")
def tmp_mesh_filename(tmp_dir: Path) -> str:
    return (tmp_dir / "shape.mesh").as_posix()


@pytest.fixture(scope="function")
def tmp_output_mesh_filename(tmp_dir: Path) -> str:
    return (tmp_dir / "shape.o.mesh").as_posix()


def test_given_rve_and_triangles_on_each_boundary_face_is_triangle_on_boundary_must_return_true() -> None:
    # Arrange

    rve = Rve(dim_x=1, dim_y=1, dim_z=1, center=(0, 0, 0))
    triangle_xmax = remesh.Triangle(np.array([0.5, -0.5, -0.5]), np.array([0.5, 0.5, -0.5]), np.array([0.5, -0.5, 0.5]),
                                    1)
    triangle_xmin = remesh.Triangle(np.array([-0.5, -0.5, -0.5]), np.array([-0.5, 0.5, -0.5]),
                                    np.array([-0.5, -0.5, 0.5]),
                                    2)
    triangle_ymax = remesh.Triangle(np.array([0.5, 0.5, -0.5]), np.array([-0.5, 0.5, -0.5]), np.array([0.5, 0.5, 0.5]),
                                    3)
    triangle_ymin = remesh.Triangle(np.array([0.5, -0.5, -0.5]), np.array([-0.5, -0.5, -0.5]),
                                    np.array([0.5, -0.5, 0.5]),
                                    4)
    triangle_zmax = remesh.Triangle(np.array([0.5, -0.5, 0.5]), np.array([0.5, 0.5, 0.5]), np.array([-0.5, -0.5, 0.5]),
                                    5)
    triangle_zmin = remesh.Triangle(np.array([0.5, -0.5, -0.5]), np.array([0.5, 0.5, -0.5]),
                                    np.array([-0.5, -0.5, -0.5]),
                                    6)

    face_triangles = [triangle_xmax, triangle_xmin, triangle_ymax, triangle_ymin, triangle_zmax, triangle_zmin]

    # Act & Assert
    assert all(remesh._is_triangle_on_boundary(triangle, rve) for triangle in face_triangles)


def test_given_rve_and_internal_triangle_is_triangle_on_boundary_must_return_false() -> None:
    # Arrange

    rve = Rve(dim_x=1, dim_y=1, dim_z=1, center=(0, 0, 0))
    internal_triangle = remesh.Triangle(np.array([0.2, -0.25, -0.25]), np.array([0.2, 0.25, -0.25]),
                               np.array([0.2, -0.25, 0.25]), 1)

    # Act & Assert
    assert not remesh._is_triangle_on_boundary(internal_triangle, rve)


@pytest.mark.parametrize(
    "shape, mesh_element_size",
    [
        (Box(), 1.0),
        (
                Tpms(
                    surface_function=tpms.gyroid,
                    type_part="sheet",
                    thickness=0.05,
                ),
                0.05,
        ),
    ],
)
def test_given_periodic_mesh_remesh_keeping_periodicity_for_fem_must_maintain_periodicity(
        shape: BasicGeometry,
        mesh_element_size: float,
        tmp_step_filename: str,
        tmp_mesh_filename: str,
        tmp_output_mesh_filename: str,
) -> None:
    # Arrange
    rve = Rve(dim_x=1, dim_y=1, dim_z=1, center=(0, 0, 0))
    cad_geometry: cq.Shape = shape.generate()
    phase = Phase(cad_geometry)

    cq.exporters.export(cad_geometry, tmp_step_filename)

    meshPeriodic(
        tmp_step_filename,
        rve,
        [phase],
        size=mesh_element_size,
        order=1,
        output_file=tmp_mesh_filename,
        mshFileVersion=4,
    )

    # Act
    remesh.remesh_keeping_periodicity_for_fem(
        tmp_mesh_filename, rve, tmp_output_mesh_filename, hgrad=1.1
    )

    # Assert
    if USE_MMG == True:
        assert _is_periodic(_get_mesh_nodes_coords(tmp_output_mesh_filename))