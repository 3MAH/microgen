"""
Mesh using gmsh
"""

from __future__ import annotations

from typing import Iterator

import gmsh
import numpy as np
import numpy.typing as npt

from .phase import Phase
from .rve import Rve

_DIM_COUNT = 3
_BOUNDS_COUNT = 2
_Point3D = np.ndarray


class OutputMeshNotPeriodicError(Exception):
    """Raised when output mesh from meshPeriodic is not periodic"""


def mesh(
    mesh_file: str,
    listPhases: list[Phase],
    size: float,
    order: int,
    output_file: str = "Mesh.msh",
    mshFileVersion: int = 4,
) -> None:
    """
    Meshes step file with gmsh with list of phases management

    :param mesh_file: step file to mesh
    :param listPhases: list of phases to mesh
    :param size: mesh size constraint (see: `gmsh.model.mesh.setSize(dimTags, size)`_)
    :param order: see `gmsh.model.mesh.setOrder(order)`_
    :param output_file: output file (.msh, .vtk)
    :param mshFileVersion: gmsh file version

    .. _gmsh.model.mesh.setOrder(order): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L1688
    .. _gmsh.model.mesh.setSize(dimTags, size): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L3140
    """
    _initialize_mesh(mesh_file, listPhases, order, mshFileVersion)
    _finalize_mesh(size, output_file)


def meshPeriodic(
    mesh_file: str,
    rve: Rve,
    listPhases: list[Phase],
    size: float,
    order: int,
    output_file: str = "MeshPeriodic.msh",
    mshFileVersion: int = 4,
    tol: float = 1e-8,
) -> None:
    """
    Meshes periodic geometries with gmsh

    :param mesh_file: step file to mesh
    :param rve: RVE for periodicity
    :param listPhases: list of phases to mesh
    :param size: mesh size constraint (see: `gmsh.model.mesh.setSize(dimTags, size)`_)
    :param order: see `gmsh.model.mesh.setOrder(order)`_
    :param output_file: output file (.msh, .vtk)
    :param mshFileVersion: gmsh file version
    :param tol: tolerance for periodicity check

    .. _gmsh.model.mesh.setOrder(order): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L1688
    .. _gmsh.model.mesh.setSize(dimTags, size): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L3140
    """
    _initialize_mesh(mesh_file, listPhases, order, mshFileVersion)
    _set_periodic(rve)
    _finalize_mesh(size, output_file)
    _check_output_mesh_periodicity(output_file, tol)


def is_periodic(
    nodes_coords: npt.NDArray[np.float_], tol: float = 1e-8, dim: int = 3
) -> bool:
    """Checks whether a mesh is periodic, given its nodes' coordinates

    :param nodes_coords: list of nodes coordinates of the analyzed mesh
    :param tol: tolerance
    :param dim: mesh dimension
    """
    # bounding box
    xmax = np.max(nodes_coords[:, 0])
    xmin = np.min(nodes_coords[:, 0])
    ymax = np.max(nodes_coords[:, 1])
    ymin = np.min(nodes_coords[:, 1])
    if dim == 3:
        zmax = np.max(nodes_coords[:, 2])
        zmin = np.min(nodes_coords[:, 2])

    # extract face nodes
    face_xm = np.where(np.abs(nodes_coords[:, 0] - xmin) < tol)[0]
    face_xp = np.where(np.abs(nodes_coords[:, 0] - xmax) < tol)[0]

    if dim > 1:
        face_ym = np.where(np.abs(nodes_coords[:, 1] - ymin) < tol)[0]
        face_yp = np.where(np.abs(nodes_coords[:, 1] - ymax) < tol)[0]

    if dim > 2:  # or dim == 3
        face_zm = np.where(np.abs(nodes_coords[:, 2] - zmin) < tol)[0]
        face_zp = np.where(np.abs(nodes_coords[:, 2] - zmax) < tol)[0]

        # sort adjacent faces to ensure node correspondence
    if nodes_coords.shape[1] == 2:  # 2D mesh
        face_xm = face_xm[np.argsort(nodes_coords[face_xm, 1])]
        face_xp = face_xp[np.argsort(nodes_coords[face_xp, 1])]
        if dim > 1:
            face_ym = face_ym[np.argsort(nodes_coords[face_ym, 0])]
            face_yp = face_yp[np.argsort(nodes_coords[face_yp, 0])]

    elif nodes_coords.shape[1] > 2:
        decimal_round = int(-np.log10(tol) - 1)

        def _sort_dim(indices: np.ndarray, dim_a: int, dim_b: int) -> np.ndarray:
            return indices[
                np.lexsort(
                    (
                        nodes_coords[indices, dim_a],
                        nodes_coords[indices, dim_b].round(decimal_round),
                    )
                )
            ]

        face_xm = _sort_dim(face_xm, dim_a=1, dim_b=2)
        face_xp = _sort_dim(face_xp, dim_a=1, dim_b=2)
        if dim > 1:
            face_ym = _sort_dim(face_ym, dim_a=0, dim_b=2)
            face_yp = _sort_dim(face_yp, dim_a=0, dim_b=2)
        if dim > 2:
            face_zm = _sort_dim(face_zm, dim_a=0, dim_b=1)
            face_zp = _sort_dim(face_zp, dim_a=0, dim_b=1)

    # ==========================
    # test if mesh is periodic:
    # ==========================

    # test if same number of nodes in adjacent faces
    if len(face_xm) != len(face_xp):
        return False
    if dim > 1 and len(face_ym) != len(face_yp):
        return False
    if dim > 2 and (len(face_zm) != len(face_zp)):
        return False

    # check nodes position
    if (nodes_coords[face_xp, 1:] - nodes_coords[face_xm, 1:] > tol).any():
        return False
    if (
        dim > 1
        and (nodes_coords[face_yp, ::2] - nodes_coords[face_ym, ::2] > tol).any()
    ):
        return False
    if dim > 2 and (nodes_coords[face_zp, :2] - nodes_coords[face_zm, :2] > tol).any():
        return False

    return True


def _generate_list_tags(list_phases: list[Phase]) -> list[list[int]]:
    list_tags: list[list[int]] = []
    start: int = 1
    for phase in list_phases:
        stop = start + len(phase.solids)
        list_tags.append(list(range(start, stop)))
        start = stop
    return list_tags


def _generate_list_dim_tags(list_phases: list[Phase]) -> list[tuple[int, int]]:
    nb_tags = sum(len(phase.solids) for phase in list_phases)
    return [(_DIM_COUNT, tag) for tag in range(1, nb_tags + 1)]


def _initialize_mesh(
    mesh_file: str,
    list_phases: list[Phase],
    order: int,
    msh_file_version: int = 4,
) -> None:
    gmsh.initialize()
    gmsh.option.setNumber(
        name="General.Verbosity", value=1
    )  # this would still print errors, but not warnings

    gmsh.model.mesh.setOrder(order=order)
    gmsh.option.setNumber(name="Mesh.MshFileVersion", value=msh_file_version)
    gmsh.model.occ.importShapes(fileName=mesh_file, highestDimOnly=True)

    list_dim_tags = _generate_list_dim_tags(list_phases)
    if len(list_dim_tags) > 1:
        gmsh.model.occ.fragment(
            objectDimTags=list_dim_tags[:-1], toolDimTags=[list_dim_tags[-1]]
        )

    gmsh.model.occ.synchronize()

    list_tags = _generate_list_tags(list_phases)
    for i, tags in enumerate(list_tags):
        ps_i = gmsh.model.addPhysicalGroup(dim=_DIM_COUNT, tags=tags)
        gmsh.model.setPhysicalName(dim=_DIM_COUNT, tag=ps_i, name="Mat" + str(i))


def _finalize_mesh(
    size: float,
    output_file: str = "Mesh.msh",
) -> None:
    list_dim_tags: list[tuple[int, int]] = gmsh.model.getEntities()
    gmsh.model.mesh.setSize(dimTags=list_dim_tags, size=size)
    gmsh.model.mesh.generate(dim=_DIM_COUNT)
    gmsh.write(fileName=output_file)
    gmsh.finalize()


def _check_output_mesh_periodicity(output_mesh_file: str, tol: float = 1e-8) -> None:
    gmsh.initialize()
    gmsh.open(output_mesh_file)
    dimension = gmsh.model.getDimension()
    _, nodes_coords, _ = gmsh.model.mesh.getNodes(dim=dimension)
    nodes_coords = nodes_coords.reshape((-1, dimension))
    check_periodicity = is_periodic(nodes_coords, tol, dimension)
    gmsh.finalize()

    if not check_periodicity:
        raise OutputMeshNotPeriodicError(
            "Something went wrong: output mesh from meshPeriodic is not periodic."
            "\n Try changing tolerance value or mesh size parameter"
        )


def _set_periodic(rve: Rve) -> None:
    for axis in range(_DIM_COUNT):
        _set_periodic_on_axis(rve, axis)


def _iter_bounding_boxes(
    minimum: _Point3D, maximum: _Point3D, eps: float
) -> Iterator[tuple[np.ndarray, int]]:
    entities: list[tuple[int, int]] = gmsh.model.getEntitiesInBoundingBox(
        *np.subtract(minimum, eps), *np.add(maximum, eps), dim=2
    )

    for dim, tag in entities:
        bounds = np.asarray(gmsh.model.getBoundingBox(dim, tag)).reshape(
            (_BOUNDS_COUNT, _DIM_COUNT)
        )
        yield bounds, tag


def _get_deltas(rve: Rve) -> np.ndarray:  # To add as a property of Rve?
    return np.array([rve.dx, rve.dy, rve.dz])


def _get_min(rve: Rve) -> np.ndarray:  # To add as a property of Rve?
    return np.array([rve.x_min, rve.y_min, rve.z_min])


def _get_max(rve: Rve) -> np.ndarray:  # To add as a property of Rve?
    return np.array([rve.x_max, rve.y_max, rve.z_max])


def _iter_matching_bounding_boxes(rve: Rve, axis: int) -> Iterator[tuple[int, int]]:
    deltas = _get_deltas(rve)
    eps: float = 1.0e-3 * min(deltas)
    minimum = _get_min(rve)
    maximum = _get_max(rve)
    maximum[axis] = minimum[axis]

    # Get all the entities on the surface m (minimum value on axis, i.e. Xm, Ym or Zm)
    for bounds_min, tag_min in _iter_bounding_boxes(minimum, maximum, eps):
        # Translate the minimal bounds into the maximum value on axis surface
        bounds_min[:, axis] += rve.dim[axis]

        # Get all the entities on the corresponding surface (i.e. Xp, Yp or Zp)
        for bounds_max, tag_max in _iter_bounding_boxes(
            bounds_min[0], bounds_min[1], eps
        ):
            if np.all(np.abs(np.subtract(bounds_max, bounds_min)) < eps):
                yield tag_min, tag_max


def _set_periodic_on_axis(rve: Rve, axis: int) -> None:
    translation_matrix = np.eye(_DIM_COUNT + 1)
    translation_matrix[axis, _DIM_COUNT] = rve.dim[axis]
    translation: list[float] = list(translation_matrix.flatten())

    for tag_min, tag_max in _iter_matching_bounding_boxes(rve, axis):
        gmsh.model.mesh.setPeriodic(
            dim=2, tags=[tag_max], tagsMaster=[tag_min], affineTransform=translation
        )
