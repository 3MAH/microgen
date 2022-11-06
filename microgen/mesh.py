"""
Mesh using gmsh
"""
from typing import Iterator

import gmsh
import numpy as np

from .phase import Phase
from .rve import Rve

_DIM_COUNT = 3
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
    _init_mesh(mesh_file, listPhases, order, mshFileVersion)
    _save_mesh(size, output_file)


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
    _init_mesh(mesh_file, listPhases, order, mshFileVersion)
    _compute_periodicity(rve)
    _save_mesh(size, output_file)


def _generate_list_tags(listPhases: list[Phase]) -> list[list[int]]:
    listTags: list[list[int]] = []
    index: int = 0
    for phase in listPhases:
        temp: list[int] = []
        for _ in phase.solids:
            index += 1
            temp.append(index)
        listTags.append(temp)
    return listTags


def _init_mesh(
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

    flatListSolids = [solid for phase in listPhases for solid in phase.solids]
    nbTags = len(flatListSolids)
    flatListTags = list(range(1, nbTags + 1, 1))

    listTags = _generate_list_tags(listPhases)

    listDimTags = [(3, tag) for tag in flatListTags]

    gmsh.model.occ.importShapes(fileName=mesh_file, highestDimOnly=True)

    if len(listDimTags) > 1:
        gmsh.model.occ.fragment(
            objectDimTags=listDimTags[:-1], toolDimTags=[listDimTags[-1]]
        )

    gmsh.model.occ.synchronize()

    for i, tag in enumerate(listTags):
        ps_i = gmsh.model.addPhysicalGroup(dim=3, tags=tag)
        gmsh.model.setPhysicalName(dim=3, tag=ps_i, name="Mat" + str(i))


def _save_mesh(
        size: float,
        output_file: str = "Mesh.msh",
) -> None:
    p = gmsh.model.getEntities()
    gmsh.model.mesh.setSize(dimTags=p, size=size)
    gmsh.model.mesh.generate(dim=3)
    gmsh.write(fileName=output_file)
    gmsh.finalize()


def _compute_periodicity(rve: Rve) -> None:
    for axis in range(3):
        _compute_periodicity_on_axis(rve, axis)


def _iter_bounding_boxes(minimum: _Point3D, maximum: _Point3D, eps: float) -> Iterator[tuple[np.ndarray, int]]:
    entities: list[tuple[int, int]] = gmsh.model.getEntitiesInBoundingBox(
        *np.subtract(minimum, eps),
        *np.add(maximum, eps),
        dim=2
    )
    for dim, tag in entities:
        bounds = np.asarray(gmsh.model.getBoundingBox(
            dim, tag
        )).reshape((2, _DIM_COUNT))
        yield bounds, tag


def _compute_periodicity_on_axis(rve: Rve, axis: int) -> None:
    translation_matrix = np.eye(4)
    translation_matrix[axis, 3] = rve.delta[axis]
    translation = list(translation_matrix.flatten())

    eps = 1.0e-3 * min(rve.delta)

    minimum = np.zeros(_DIM_COUNT)
    maximum = np.array(rve.delta)
    maximum[axis] = 0.

    for bounds_min, tag_min in _iter_bounding_boxes(minimum, maximum, eps):
        bounds_min[:, axis] += 1
        for bounds_max, tag_max in _iter_bounding_boxes(bounds_min[0], bounds_min[1], eps):
            bounds_max[:, axis] -= 1
            if (
                    np.all(np.abs(np.subtract(bounds_max, bounds_min)) < eps)
            ):
                gmsh.model.mesh.setPeriodic(
                    dim=2,
                    tags=[tag_max],
                    tagsMaster=[tag_min],
                    affineTransform=translation
                )
