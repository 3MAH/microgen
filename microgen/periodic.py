"""Periodic function to cut a shape periodically according to a RVE."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import cadquery as cq

from .operations import fuseShapes
from .phase import Phase

if TYPE_CHECKING:
    from .rve import Rve

_NO_INTERSECTION = 0
_FACE = 1
_EDGE = 2
_CORNER = 3

_FACES = ["x-", "x+", "y-", "y+", "z-", "z+"]
_DIRECTION = {
    "x-": (-1, 0, 0),
    "x+": (1, 0, 0),
    "y-": (0, -1, 0),
    "y+": (0, 1, 0),
    "z-": (0, 0, -1),
    "z+": (0, 0, 1),
}
FACE_DIR = {"x-": ">X", "x+": "<X", "y-": ">Y", "y+": "<Y", "z-": ">Z", "z+": "<Z"}
INV_FACE_DIR = {
    "x-": "<X",
    "x+": ">X",
    "y-": "<Y",
    "y+": ">Y",
    "z-": "<Z",
    "z+": ">Z",
}


def _detect_intersected_faces(
    wk_plane: cq.Workplane,
    rve: Rve,
) -> tuple[list[str], dict[str, cq.Face], dict[str, cq.Workplane]]:
    base_pnt = {
        "x-": (rve.min_point[0], 0, 0),
        "x+": (rve.max_point[0], 0, 0),
        "y-": (0, rve.min_point[1], 0),
        "y+": (0, rve.max_point[1], 0),
        "z-": (0, 0, rve.min_point[2]),
        "z+": (0, 0, rve.max_point[2]),
    }
    intersected_faces: list[str] = []
    rve_planes: dict[str, cq.Face] = {}
    partitions: dict[str, cq.Workplane] = {}
    for face in _FACES:
        rve_planes[face] = cq.Face.makePlane(
            basePnt=base_pnt[face],
            dir=_DIRECTION[face],
        )
        partitions[face] = wk_plane.split(cq.Workplane().add(rve_planes[face]))

        # detection of intersected faces
        if len(partitions[face].solids().all()) > 1:
            intersected_faces.append(face)

    # check if the object is intersecting two opposite faces
    for axis in "xyz":
        face_p = f"{axis}+"
        face_m = f"{axis}-"
        if face_p in intersected_faces and face_m in intersected_faces:
            intersected_faces.remove(face_p)
            intersected_faces.remove(face_m)
            warnings.warn(
                f"Object intersecting {face_p} and {face_m} faces: "
                "not doing anything in this direction",
                stacklevel=2,
            )
    return intersected_faces, rve_planes, partitions


def _intersect_face(
    face: str,
    rve: Rve,
    partitions: dict[str, cq.Workplane],
) -> list[cq.Workplane]:
    translate = {
        "x-": (rve.dim[0], 0, 0),
        "x+": (-rve.dim[0], 0, 0),
        "y-": (0, rve.dim[1], 0),
        "y+": (0, -rve.dim[1], 0),
        "z-": (0, 0, rve.dim[2]),
        "z+": (0, 0, -rve.dim[2]),
    }
    inside = partitions[face].solids(FACE_DIR[face]).intersect(rve.box)
    outside = (
        partitions[face]
        .solids(INV_FACE_DIR[face])
        .translate(translate[face])
        .intersect(rve.box)
    )
    return [inside, outside]


def _intersect_edge(
    faces: tuple[str, str],
    rve: Rve,
    partitions: dict[str, cq.Workplane],
    rve_planes: dict[str, cq.Face],
) -> list[cq.Workplane]:
    translate = {
        "x-": (rve.dim[0], 0, 0),
        "x+": (-rve.dim[0], 0, 0),
        "y-": (0, rve.dim[1], 0),
        "y+": (0, -rve.dim[1], 0),
        "z-": (0, 0, rve.dim[2]),
        "z+": (0, 0, -rve.dim[2]),
    }
    f_0, f_1 = faces
    part = partitions[f_0].solids(FACE_DIR[f_0]).split(rve_planes[f_1])
    periodic_object = [
        part.solids(FACE_DIR[f_1]).intersect(rve.box),
        part.solids(INV_FACE_DIR[f_1]).translate(translate[f_1]).intersect(rve.box),
    ]
    part = partitions[f_0].solids(INV_FACE_DIR[f_0]).split(rve_planes[f_1])
    periodic_object.append(
        part.solids(FACE_DIR[f_1]).translate(translate[f_0]).intersect(rve.box),
    )
    tslt = (
        translate[f_0][0] + translate[f_1][0],
        translate[f_0][1] + translate[f_1][1],
        translate[f_0][2] + translate[f_1][2],
    )
    periodic_object.append(
        part.solids(INV_FACE_DIR[f_1]).translate(tslt).intersect(rve.box),
    )
    return periodic_object


def _intersect_corner(
    faces: tuple[str, str, str],
    rve: Rve,
    partitions: dict[str, cq.Workplane],
    rve_planes: dict[str, cq.Face],
) -> list[cq.Workplane]:
    translate = {
        "x-": (rve.dim[0], 0, 0),
        "x+": (-rve.dim[0], 0, 0),
        "y-": (0, rve.dim[1], 0),
        "y+": (0, -rve.dim[1], 0),
        "z-": (0, 0, rve.dim[2]),
        "z+": (0, 0, -rve.dim[2]),
    }
    f_0, f_1, f_2 = faces

    new_part = (
        partitions[f_0]
        .solids(FACE_DIR[f_0])
        .split(rve_planes[f_1])
        .solids(FACE_DIR[f_1])
    )
    periodic_object = [
        new_part.solids(FACE_DIR[f_2]).intersect(rve.box),
        new_part.solids(INV_FACE_DIR[f_2]).translate(translate[f_2]).intersect(rve.box),
    ]
    new_part = (
        partitions[f_0]
        .solids(FACE_DIR[f_0])
        .split(rve_planes[f_1])
        .solids(INV_FACE_DIR[f_1])
    )
    periodic_object.extend(
        (
            new_part.solids(FACE_DIR[f_1]).translate(translate[f_1]).intersect(rve.box),
            new_part.solids(INV_FACE_DIR[f_1])
            .translate((0, translate[f_1][1], translate[f_2][2]))
            .intersect(rve.box),
        ),
    )
    new_part = (
        partitions[f_0]
        .solids(INV_FACE_DIR[f_0])
        .split(rve_planes[f_1])
        .solids(FACE_DIR[f_1])
    )
    periodic_object.extend(
        (
            new_part.solids(FACE_DIR[f_2]).translate(translate[f_0]).intersect(rve.box),
            new_part.solids(INV_FACE_DIR[f_2])
            .translate((translate[f_0][0], 0, translate[f_2][2]))
            .intersect(rve.box),
        ),
    )
    new_part = (
        partitions[f_0]
        .solids(INV_FACE_DIR[f_0])
        .split(rve_planes[f_1])
        .solids(INV_FACE_DIR[f_1])
    )
    periodic_object.extend(
        (
            new_part.solids(FACE_DIR[f_2])
            .translate((translate[f_0][0], translate[f_1][1], 0))
            .intersect(rve.box),
            new_part.solids(INV_FACE_DIR[f_2])
            .translate((translate[f_0][0], translate[f_1][1], translate[f_2][2]))
            .intersect(rve.box),
        ),
    )
    return periodic_object


def periodic_split_and_translate(phase: Phase, rve: Rve) -> Phase:
    """Rearrange phase periodically according to the rve.

    :param phase: Phase to cut periodically
    :param rve: RVE for periodicity

    :return phase: resulting phase
    """
    wk_plane = cq.Workplane().add(phase.solids)  # shape to cut

    intersected_faces, rve_planes, partitions = _detect_intersected_faces(wk_plane, rve)

    periodic_object: list[cq.Workplane] = []
    if len(intersected_faces) == _NO_INTERSECTION:
        periodic_object.append(wk_plane)

    elif len(intersected_faces) == _FACE:
        periodic_object.extend(
            _intersect_face(
                face=intersected_faces[0],
                rve=rve,
                partitions=partitions,
            ),
        )

    elif len(intersected_faces) == _EDGE:
        periodic_object.extend(
            _intersect_edge(
                faces=(intersected_faces[0], intersected_faces[1]),
                rve=rve,
                partitions=partitions,
                rve_planes=rve_planes,
            ),
        )

    elif len(intersected_faces) == _CORNER:
        periodic_object.extend(
            _intersect_corner(
                faces=(
                    intersected_faces[0],
                    intersected_faces[1],
                    intersected_faces[2],
                ),
                rve=rve,
                partitions=partitions,
                rve_planes=rve_planes,
            ),
        )

    list_solids = [wp.val().Solids() for wp in periodic_object]
    flat_list = [solid.copy() for solids in list_solids for solid in solids]
    to_fuse = [cq.Shape(solid.wrapped) for solid in flat_list]
    shape = fuseShapes(cqShapeList=to_fuse, retain_edges=False)

    return Phase(shape=shape)


def periodic(phase: Phase, rve: Rve) -> Phase:
    """See periodic_split_and_translate.

    Deprecated in favor of periodic_split_and_translate.
    """
    warnings.warn(
        "periodic is deprecated, use periodic_split_and_translate instead",
        DeprecationWarning,
        stacklevel=2,
    )
    return periodic_split_and_translate(phase, rve)
