import subprocess
import pytest
from microgen import Rve, remesh, Tpms, BoxMesh
from microgen.shape.surface_functions import gyroid
import numpy as np
import numpy.typing as npt
import gmsh
import meshio
import pyvista as pv
from pathlib import Path

_MESH_DIM = 3
_BOUNDARY_DIM = 2

USE_MMG = False
try:
    subprocess.check_output("mmg3d_O3", stderr=subprocess.STDOUT)
    USE_MMG= True    
except(subprocess.CalledProcessError, FileNotFoundError):
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
def tmp_mesh_filename(tmp_dir: Path) -> str:
    return (tmp_dir / "shape.mesh").as_posix()

@pytest.fixture(scope="function")
def tmp_vtk_filename(tmp_dir: Path) -> str:
    return (tmp_dir / "shape.vtk").as_posix()


@pytest.fixture(scope="function")
def tmp_output_mesh_filename(tmp_dir: Path) -> str:
    return (tmp_dir / "shape.o.mesh").as_posix()

@pytest.fixture(scope='session')
def rve() -> Rve:
    return Rve(dim_x=1.0, dim_y=1.0, dim_z=1.0, center=(0.5, 0.5, 0.5))

@pytest.fixture(scope="function")
def box(rve: Rve) -> BoxMesh:
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

    elements_dict = {pv.CellType.TETRA: np.array([[11, 6, 0, 14],
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

    mesh = BoxMesh(nodes_array, elements_dict)
    mesh.rve = rve

    return mesh

@pytest.fixture(scope='function')
def gyroid(rve: Rve) -> BoxMesh:
    gyroid = Tpms(surface_function=gyroid, offset=1.0).generateVtk(type_part="sheet")
    gyroid_boxmesh = BoxMesh.from_pyvista(gyroid)
    return gyroid_boxmesh



@pytest.mark.parametrize(
    "shape",
    [
        "box",
        "gyroid",
    ]
)
def test_given_periodic_mesh_remesh_keeping_periodicity_for_fem_must_maintain_periodicity(
        shape: BoxMesh,
        request,
        rve: Rve,
        tmp_mesh_filename: str,
        tmp_vtk_filename: str,
        tmp_output_mesh_filename: str,
) -> None:
    if USE_MMG == True:
        # Arrange

        request.getfixturevalue(shape).to_pyvista().save(tmp_vtk_filename)
        initial_mesh = meshio.read(tmp_vtk_filename)
        initial_mesh.write(tmp_mesh_filename)

        # Act
        
        remesh.remesh_keeping_periodicity_for_fem(
            shape, tmp_output_mesh_filename, hgrad=1.05
        )

        # Assert
        assert _is_periodic(_get_mesh_nodes_coords(tmp_mesh_filename)) and _is_periodic(_get_mesh_nodes_coords(tmp_output_mesh_filename))