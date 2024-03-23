"""
Periodic function to cut a shape periodically according to a RVE
"""

import warnings
from typing import Dict, List, Tuple

import cadquery as cq
import numpy as np

from .operations import fuseShapes
from .phase import Phase
from .rve import Rve

FACES = ["x-", "x+", "y-", "y+", "z-", "z+"]

DIRECTION = {
    "x-": (-1, 0, 0),
    "x+": (1, 0, 0),
    "y-": (0, -1, 0),
    "y+": (0, 1, 0),
    "z-": (0, 0, -1),
    "z+": (0, 0, 1),
}

FACE_DIR = {
    "x-": ">X",
    "x+": "<X",
    "y-": ">Y",
    "y+": "<Y",
    "z-": ">Z",
    "z+": "<Z",
}
INVERSE_FACE_DIR = {
    "x-": "<X",
    "x+": ">X",
    "y-": "<Y",
    "y+": ">Y",
    "z-": "<Z",
    "z+": ">Z",
}


def periodic(phase: Phase, rve: Rve) -> Phase:
    """
    Rearrange phase periodically according to the rve

    :param phase: Phase to cut periodically
    :param rve: RVE for periodicity

    :return phase: resulting phase
    """

    periodic_object: List[cq.Workplane]
    wk_plane = cq.Workplane().add(phase.solids)  # shape to cut

    basePnt = {
        "x-": (rve.min_point[0], 0, 0),
        "x+": (rve.max_point[0], 0, 0),
        "y-": (0, rve.min_point[1], 0),
        "y+": (0, rve.max_point[1], 0),
        "z-": (0, 0, rve.min_point[2]),
        "z+": (0, 0, rve.max_point[2]),
    }

    translate = {
        "x-": (rve.dim[0], 0, 0),
        "x+": (-rve.dim[0], 0, 0),
        "y-": (0, rve.dim[1], 0),
        "y+": (0, -rve.dim[1], 0),
        "z-": (0, 0, rve.dim[2]),
        "z+": (0, 0, -rve.dim[2]),
    }

    intersected_faces, planes, partitions = _get_periodic_splitting_infos(
        base_pnt=basePnt, wk_plane=wk_plane
    )

    intersected_faces = _remove_opposite_intersected_faces(intersected_faces)

    if not intersected_faces:  # if no intersected faces = nothing to do
        periodic_object = [wk_plane]

    elif len(intersected_faces) == 1:  # one face intersected
        periodic_object = _intersect_face(
            intersected_faces=intersected_faces,
            partitions=partitions,
            rve=rve,
            translate=translate,
        )
    elif len(intersected_faces) == 2:  # two faces intersected (edge)
        periodic_object = _intersect_edge(
            intersected_faces=intersected_faces,
            partitions=partitions,
            rve=rve,
            translate=translate,
            planes=planes,
        )
    elif len(intersected_faces) == 3:  # three faces intersected (corner)
        periodic_object = _intersect_corner(
            intersected_faces=intersected_faces,
            partitions=partitions,
            rve=rve,
            translate=translate,
            planes=planes,
        )
    solids = [solid.copy() for wp in periodic_object for solid in wp.val().Solids()]
    shape = fuseShapes(
        cqShapeList=[cq.Shape(solid.wrapped) for solid in solids],
        retain_edges=False,
    )

    return Phase(shape=shape)


def _get_periodic_splitting_infos(
    base_pnt: Dict[str, np.ndarray],
    wk_plane: cq.Workplane,
) -> Tuple[
    List[str],
    Dict[str, cq.Face],
    Dict[str, cq.Workplane],
]:
    intersected_faces: List[str] = []
    planes: Dict[str, cq.Face] = {}
    partitions: Dict[str, cq.Workplane] = {}
    for face in FACES:
        planes[face] = cq.Face.makePlane(
            basePnt=base_pnt[face], dir=DIRECTION[face]
        )  # planes composing the Rve box
        partitions[face] = wk_plane.split(
            cq.Workplane().add(planes[face])
        )  # each partition cuts the object by the corresponding plane

        # detection of intersected faces
        if len(partitions[face].solids().all()) > 1:
            intersected_faces.append(face)

    return intersected_faces, planes, partitions


def _remove_opposite_intersected_faces(intersected_faces: List[str]) -> List[str]:
    for axis in ["x", "y", "z"]:
        if f"{axis}-" in intersected_faces and f"{axis}+" in intersected_faces:
            intersected_faces.remove(f"{axis}-")
            intersected_faces.remove(f"{axis}+")
            warnings.warn(
                f"Object intersecting {axis}+ and {axis}- faces: \
                    not doing anything in this direction"
            )
    return intersected_faces


def _intersect_face(
    intersected_faces: List[str],
    partitions: Dict[str, cq.Workplane],
    rve: Rve,
    translate: Dict[str, Tuple[float, float, float]],
) -> List[cq.Workplane]:
    f_0 = intersected_faces[0]
    return [
        partitions[f_0].solids(FACE_DIR[f_0]).intersect(rve.box),  # add the part of the
        partitions[f_0]
        .solids(INVERSE_FACE_DIR[f_0])  # translate the outside part of
        .translate(translate[f_0])  # the object in the rve and add
        .intersect(rve.box),
    ]


def _intersect_edge(
    intersected_faces: List[str],
    partitions: Dict[str, cq.Workplane],
    rve: Rve,
    translate: Dict[str, Tuple[float, float, float]],
    planes: Dict[str, cq.Face],
) -> List[cq.Workplane]:
    f_0 = intersected_faces[0]
    f_1 = intersected_faces[1]

    part = (
        cq.Workplane()
        .add(partitions[f_0].solids(FACE_DIR[f_0]))
        .split(cq.Workplane().add(planes[f_1]))
    )
    periodic_object = [
        part.solids(FACE_DIR[f_1]).intersect(rve.box),
        part.solids(INVERSE_FACE_DIR[f_1]).translate(translate[f_1]).intersect(rve.box),
    ]
    tslt = tuple(translate[f_0][i] + translate[f_1][i] for i in range(3))
    part = (
        cq.Workplane()
        .add(partitions[f_0].solids(INVERSE_FACE_DIR[f_0]))
        .split(cq.Workplane().add(planes[f_1]))
    )
    periodic_object.extend(
        (
            part.solids(FACE_DIR[f_1]).translate(translate[f_0]).intersect(rve.box),
            part.solids(INVERSE_FACE_DIR[f_1]).translate(tslt).intersect(rve.box),
        )
    )
    return periodic_object


def _intersect_corner(
    intersected_faces: List[str],
    partitions: Dict[str, cq.Workplane],
    rve: Rve,
    translate: Dict[str, Tuple[float, float, float]],
    planes: Dict[str, cq.Face],
) -> List[cq.Workplane]:
    f_0 = intersected_faces[0]
    f_1 = intersected_faces[1]
    f_2 = intersected_faces[2]

    new_part = (
        cq.Workplane()
        .add(partitions[f_0].solids(FACE_DIR[f_0]))
        .split(cq.Workplane().add(planes[f_1]))
        .solids(FACE_DIR[f_1])
    )
    periodic_object = [
        new_part.solids(FACE_DIR[f_2]).intersect(rve.box),
        new_part.solids(INVERSE_FACE_DIR[f_2])
        .translate(translate[f_2])
        .intersect(rve.box),
    ]
    new_part = (
        cq.Workplane()
        .add(partitions[f_0].solids(FACE_DIR[f_0]))
        .split(cq.Workplane().add(planes[f_1]))
        .solids(INVERSE_FACE_DIR[f_1])
    )
    periodic_object.extend(
        (
            new_part.solids(FACE_DIR[f_1]).translate(translate[f_1]).intersect(rve.box),
            new_part.solids(INVERSE_FACE_DIR[f_1])
            .translate((0, translate[f_1][1], translate[f_2][2]))
            .intersect(rve.box),
        )
    )
    new_part = (
        cq.Workplane()
        .add(partitions[f_0].solids(INVERSE_FACE_DIR[f_0]))
        .split(cq.Workplane().add(planes[f_1]))
        .solids(FACE_DIR[f_1])
    )
    periodic_object.extend(
        (
            new_part.solids(FACE_DIR[f_2]).translate(translate[f_0]).intersect(rve.box),
            new_part.solids(INVERSE_FACE_DIR[f_2])
            .translate((translate[f_0][0], 0, translate[f_2][2]))
            .intersect(rve.box),
        )
    )
    new_part = (
        cq.Workplane()
        .add(partitions[f_0].solids(INVERSE_FACE_DIR[f_0]))
        .split(cq.Workplane().add(planes[f_1]))
        .solids(INVERSE_FACE_DIR[f_1])
    )
    periodic_object.extend(
        (
            new_part.solids(FACE_DIR[f_2])
            .translate((translate[f_0][0], translate[f_1][1], 0))
            .intersect(rve.box),
            new_part.solids(INVERSE_FACE_DIR[f_2])
            .translate((translate[f_0][0], translate[f_1][1], translate[f_2][2]))
            .intersect(rve.box),
        )
    )
    return periodic_object
