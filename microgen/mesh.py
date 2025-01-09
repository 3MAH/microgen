"""Mesh using gmsh."""

from __future__ import annotations

import itertools
import warnings
from typing import TYPE_CHECKING, Iterator

import gmsh
import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from .phase import Phase
    from .rve import Rve

_DIM_COUNT = 3
_BOUNDS_COUNT = 2
_Point3D = npt.NDArray[np.float64]
_DIM_2D = 2
_DIM_3D = 3


class OutputMeshNotPeriodicError(Exception):
    """Raised when output mesh from meshPeriodic is not periodic."""


def mesh(
    mesh_file: str,
    list_phases: list[Phase] | None = None,
    size: float | None = None,
    order: int | None = None,
    output_file: str = "Mesh.msh",
    msh_file_version: int = 4,
    listPhases: list[Phase] | None = None,
    mshFileVersion: int | None = None,
) -> None:
    """Mesh step file with gmsh with list of phases management.

    :param mesh_file: step file to mesh
    :param listPhases: list of phases to mesh
    :param size: mesh size constraint (see: `gmsh.model.mesh.setSize(dimTags, size)`_)
    :param order: see `gmsh.model.mesh.setOrder(order)`_
    :param output_file: output file (.msh, .vtk)
    :param mshFileVersion: gmsh file version

    .. _gmsh.model.mesh.setOrder(order): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L1688
    .. _gmsh.model.mesh.setSize(dimTags, size): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L3140
    """
    if listPhases is not None:
        warnings.warn(
            "listPhases is deprecated, use list_phases instead",
            DeprecationWarning,
            stacklevel=2,
        )
        list_phases = listPhases

    if mshFileVersion is not None:
        warnings.warn(
            "mshFileVersion is deprecated, use msh_file_version instead",
            DeprecationWarning,
            stacklevel=2,
        )
        msh_file_version = mshFileVersion

    if size is None:
        err_msg = "size parameter must be provided"
        raise ValueError(err_msg)

    if order is None:
        err_msg = "order parameter must be provided"
        raise ValueError(err_msg)

    _initialize_mesh(mesh_file, list_phases, order, msh_file_version)
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
    """Mesh periodic geometries with gmsh.

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
    warnings.warn(
        "meshPeriodic is deprecated, use mesh_periodic instead",
        DeprecationWarning,
        stacklevel=2,
    )
    mesh_periodic(
        mesh_file,
        rve,
        listPhases,
        size,
        order,
        output_file,
        mshFileVersion,
        tol,
    )


def mesh_periodic(
    mesh_file: str,
    rve: Rve,
    list_phases: list[Phase],
    size: float,
    order: int,
    output_file: str = "MeshPeriodic.msh",
    msh_file_version: int = 4,
    tol: float = 1e-8,
) -> None:
    """Mesh periodic geometries with gmsh.

    :param mesh_file: step file to mesh
    :param rve: RVE for periodicity
    :param list_phases: list of phases to mesh
    :param size: mesh size constraint (see: `gmsh.model.mesh.setSize(dimTags, size)`_)
    :param order: see `gmsh.model.mesh.setOrder(order)`_
    :param output_file: output file (.msh, .vtk)
    :param msh_file_version: gmsh file version
    :param tol: tolerance for periodicity check

    .. _gmsh.model.mesh.setOrder(order): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L1688
    .. _gmsh.model.mesh.setSize(dimTags, size): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L3140
    """
    _initialize_mesh(mesh_file, list_phases, order, msh_file_version)
    _set_periodic(rve)
    _finalize_mesh(size, output_file)
    _check_output_mesh_periodicity(output_file, tol)


def _get_bounding_box(
    nodes_coords: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Get bounding box of a mesh.

    :param nodes_coords: list of nodes coordinates of the analyzed mesh
    """
    min_point = np.min(nodes_coords, axis=0)
    max_point = np.max(nodes_coords, axis=0)
    return min_point, max_point


def _extract_face_nodes(
    nodes_coords: npt.NDArray[np.float64],
    min_point: npt.NDArray[np.float64],
    max_point: npt.NDArray[np.float64],
    tol: float,
    dim: int,
) -> dict[str, npt.NDArray[np.int64]]:
    """Extract face nodes of a mesh.

    :param nodes_coords: list of nodes coordinates of the analyzed mesh
    :param min_point: minimum point of the bounding box
    :param max_point: maximum point of the bounding box
    :param tol: tolerance
    :param dim: mesh dimension

    :return: dictionary of face nodes as
        {
            "x-": face_xm,
            "x+": face_xp,
            "y-": face_ym,
            "y+": face_yp,
            "z-": face_zm,
            "z+": face_zp,
        }
    """
    axes = "xyz"[:dim]
    faces = {}
    for i, axis in enumerate(axes):
        for sign, point in zip("-+", [min_point, max_point]):
            face = f"{axis}{sign}"
            faces[face] = np.where(np.abs(nodes_coords[:, i] - point[i]) < tol)[0]
    return faces


def _sort_adjacent_faces_2d(
    faces: dict[str, npt.NDArray[np.int64]],
    nodes_coords: npt.NDArray[np.float64],
) -> dict[str, npt.NDArray[np.int64]]:
    """Sort adjacent faces to ensure node correspondence."""
    complementary_indices = {
        "x": 1,
        "y": 0,
    }
    for axis, sign in itertools.product("xy", "-+"):
        face = f"{axis}{sign}"
        idx = complementary_indices[axis]
        faces[face] = faces[face][np.argsort(nodes_coords[faces[face], idx])]
    return faces


def _sort_adjacent_faces_3d(
    faces: dict[str, npt.NDArray[np.int64]],
    nodes_coords: npt.NDArray[np.float64],
    tol: float,
) -> dict[str, npt.NDArray[np.int64]]:
    """Sort adjacent faces to ensure node correspondence."""
    decimal_round = int(-np.log10(tol) - 1)

    def _sort_dim(
        indices: npt.NDArray[np.int64],
        dim_a: int,
        dim_b: int,
    ) -> npt.NDArray[np.int64]:
        return indices[
            np.lexsort(
                (
                    nodes_coords[indices, dim_a],
                    nodes_coords[indices, dim_b].round(decimal_round),
                ),
            )
        ]

    complementary_indices = {
        "x": (1, 2),
        "y": (0, 2),
        "z": (0, 1),
    }
    for axis, sign in itertools.product("xyz", "-+"):
        face = f"{axis}{sign}"
        slc = complementary_indices[axis]
        faces[face] = _sort_dim(faces[face], slc[0], slc[1])
    return faces


def is_periodic(
    nodes_coords: npt.NDArray[np.float64],
    tol: float = 1e-8,
    dim: int | None = None,
) -> bool:
    """Check whether a mesh is periodic, given its nodes' coordinates.

    :param nodes_coords: list of nodes coordinates of the analyzed mesh
    :param tol: tolerance
    :param dim: mesh dimension (deprecated: unnecessary parameter)
    """
    if dim is not None:
        warnings.warn(
            "dim is deprecated, it is now inferred from the nodes' coordinates.",
            DeprecationWarning,
            stacklevel=2,
        )
    dim = nodes_coords.shape[1]
    axes = "xyz"[:dim]
    min_point, max_point = _get_bounding_box(nodes_coords)

    faces = _extract_face_nodes(nodes_coords, min_point, max_point, tol, dim)

    if dim == _DIM_2D:
        faces = _sort_adjacent_faces_2d(faces, nodes_coords)
    elif dim == _DIM_3D:
        faces = _sort_adjacent_faces_3d(faces, nodes_coords, tol)

    # test if mesh is periodic
    complementary_indices = {
        "x": slice(1, None),
        "y": slice(None, None, 2),
        "z": slice(None, 2),
    }
    for axis in axes:
        face_m = faces[f"{axis}-"]
        face_p = faces[f"{axis}+"]

        # test if same number of nodes in adjacent faces
        if len(face_m) != len(face_p):
            return False

        # check nodes position
        slc = complementary_indices[axis]
        diff = nodes_coords[face_p, slc] - nodes_coords[face_m, slc]
        if (diff > tol).any():
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
        name="General.Verbosity",
        value=1,
    )  # this would still print errors, but not warnings

    gmsh.model.mesh.setOrder(order=order)
    gmsh.option.setNumber(name="Mesh.MshFileVersion", value=msh_file_version)
    gmsh.model.occ.importShapes(fileName=mesh_file, highestDimOnly=True)

    list_dim_tags = _generate_list_dim_tags(list_phases)
    if len(list_dim_tags) > 1:
        gmsh.model.occ.fragment(
            objectDimTags=list_dim_tags[:-1],
            toolDimTags=[list_dim_tags[-1]],
        )

    gmsh.model.occ.synchronize()

    list_tags = _generate_list_tags(list_phases)
    for i, tags in enumerate(list_tags):
        ps_i = gmsh.model.addPhysicalGroup(dim=_DIM_COUNT, tags=tags)
        gmsh.model.setPhysicalName(dim=_DIM_COUNT, tag=ps_i, name=f"Mat{i}")


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
        err_msg = (
            "Something went wrong: output mesh from meshPeriodic is not periodic."
            "\n Try changing tolerance value or mesh size parameter"
        )
        raise OutputMeshNotPeriodicError(err_msg)


def _set_periodic(rve: Rve) -> None:
    for axis in range(_DIM_COUNT):
        _set_periodic_on_axis(rve, axis)


def _iter_bounding_boxes(
    minimum: _Point3D,
    maximum: _Point3D,
    eps: float,
) -> Iterator[tuple[np.ndarray, int]]:
    entities: list[tuple[int, int]] = gmsh.model.getEntitiesInBoundingBox(
        *np.subtract(minimum, eps),
        *np.add(maximum, eps),
        dim=2,
    )

    for dim, tag in entities:
        bounds = np.asarray(gmsh.model.getBoundingBox(dim, tag)).reshape(
            (_BOUNDS_COUNT, _DIM_COUNT),
        )
        yield bounds, tag


def _iter_matching_bounding_boxes(rve: Rve, axis: int) -> Iterator[tuple[int, int]]:
    eps: float = 1.0e-3 * min(rve.dim)
    minimum = rve.min_point.copy()
    maximum = rve.max_point.copy()
    maximum[axis] = minimum[axis]

    # Get all the entities on the surface m (minimum value on axis, i.e. Xm, Ym or Zm)
    for bounds_min, tag_min in _iter_bounding_boxes(minimum, maximum, eps):
        # Translate the minimal bounds into the maximum value on axis surface
        bounds_min[:, axis] += rve.dim[axis]

        # Get all the entities on the corresponding surface (i.e. Xp, Yp or Zp)
        for bounds_max, tag_max in _iter_bounding_boxes(
            bounds_min[0],
            bounds_min[1],
            eps,
        ):
            if np.all(np.abs(np.subtract(bounds_max, bounds_min)) < eps):
                yield tag_min, tag_max


def _set_periodic_on_axis(rve: Rve, axis: int) -> None:
    translation_matrix = np.eye(_DIM_COUNT + 1)
    translation_matrix[axis, _DIM_COUNT] = rve.dim[axis]
    translation: list[float] = list(translation_matrix.flatten())

    for tag_min, tag_max in _iter_matching_bounding_boxes(rve, axis):
        gmsh.model.mesh.setPeriodic(
            dim=2,
            tags=[tag_max],
            tagsMaster=[tag_min],
            affineTransform=translation,
        )
