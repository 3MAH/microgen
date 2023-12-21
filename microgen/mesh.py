"""
Mesh using gmsh
"""
from typing import Iterator

import gmsh
import numpy as np

from .phase import Phase
from .rve import Rve

_DIM_COUNT = 3
_BOUNDS_COUNT = 2
_Point3D = np.ndarray


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

    .. _gmsh.model.mesh.setOrder(order): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L1688
    .. _gmsh.model.mesh.setSize(dimTags, size): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L3140
    """
    _initialize_mesh(mesh_file, listPhases, order, mshFileVersion)
    _set_periodic(rve)
    _finalize_mesh(size, output_file)


def _generate_list_tags(listPhases: list[Phase]) -> list[list[int]]:
    listTags: list[list[int]] = []
    start: int = 1
    for phase in listPhases:
        stop = start + len(phase.solids)
        listTags.append(list(range(start, stop)))
        start = stop
    return listTags


def _generate_list_dim_tags(listPhases: list[Phase]) -> list[tuple[int, int]]:
    nbTags = sum(len(phase.solids) for phase in listPhases)
    return [(_DIM_COUNT, tag) for tag in range(1, nbTags + 1)]


def _initialize_mesh(
    mesh_file: str,
    listPhases: list[Phase],
    order: int,
    mshFileVersion: int = 4,
) -> None:
    gmsh.initialize()
    gmsh.option.setNumber(
        name="General.Verbosity", value=1
    )  # this would still print errors, but not warnings

    gmsh.model.mesh.setOrder(order=order)
    gmsh.option.setNumber(name="Mesh.MshFileVersion", value=mshFileVersion)
    gmsh.model.occ.importShapes(fileName=mesh_file, highestDimOnly=True)

    listDimTags = _generate_list_dim_tags(listPhases)
    if len(listDimTags) > 1:
        gmsh.model.occ.fragment(
            objectDimTags=listDimTags[:-1], toolDimTags=[listDimTags[-1]]
        )

    gmsh.model.occ.synchronize()

    listTags = _generate_list_tags(listPhases)
    for i, tags in enumerate(listTags):
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
        bounds_min[:, axis] += _get_deltas(rve)[axis]

        # Get all the entities on the corresponding surface (i.e. Xp, Yp or Zp)
        for bounds_max, tag_max in _iter_bounding_boxes(
            bounds_min[0], bounds_min[1], eps
        ):
            if np.all(np.abs(np.subtract(bounds_max, bounds_min)) < eps):
                yield tag_min, tag_max


def _set_periodic_on_axis(rve: Rve, axis: int) -> None:
    deltas = _get_deltas(rve)
    translation_matrix = np.eye(_DIM_COUNT + 1)
    translation_matrix[axis, _DIM_COUNT] = deltas[axis]
    translation: list[float] = list(translation_matrix.flatten())

    for tag_min, tag_max in _iter_matching_bounding_boxes(rve, axis):
        gmsh.model.mesh.setPeriodic(
            dim=2, tags=[tag_max], tagsMaster=[tag_min], affineTransform=translation
        )
