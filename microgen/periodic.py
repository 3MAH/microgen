"""
Periodic function to cut a shape periodically according to a RVE
"""

import warnings

import cadquery as cq

from .operations import fuseShapes
from .phase import Phase
from .rve import Rve


def periodic(phase: Phase, rve: Rve) -> Phase:
    """
    Rearrange phase periodically according to the rve

    :param phase: Phase to cut periodically
    :param rve: RVE for periodicity

    :return phase: resulting phase
    """

    wk_plane = cq.Workplane().add(phase.solids)  # shape to cut
    periodic_object = []  # type: list[cq.Workplane]

    faces = ["x-", "x+", "y-", "y+", "z-", "z+"]

    direction = {
        "x-": (-1, 0, 0),
        "x+": (1, 0, 0),
        "y-": (0, -1, 0),
        "y+": (0, 1, 0),
        "z-": (0, 0, -1),
        "z+": (0, 0, 1),
    }
    basePnt = {
        "x-": (rve.x_min, 0, 0),
        "x+": (rve.x_max, 0, 0),
        "y-": (0, rve.y_min, 0),
        "y+": (0, rve.y_max, 0),
        "z-": (0, 0, rve.z_min),
        "z+": (0, 0, rve.z_max),
    }

    face_dir = {"x-": ">X", "x+": "<X", "y-": ">Y", "y+": "<Y", "z-": ">Z", "z+": "<Z"}
    inverse_face_dir = {
        "x-": "<X",
        "x+": ">X",
        "y-": "<Y",
        "y+": ">Y",
        "z-": "<Z",
        "z+": ">Z",
    }

    translate = {
        "x-": (rve.dx, 0, 0),
        "x+": (-rve.dx, 0, 0),
        "y-": (0, rve.dy, 0),
        "y+": (0, -rve.dy, 0),
        "z-": (0, 0, rve.dz),
        "z+": (0, 0, -rve.dz),
    }

    intersected_faces = []  # type: list[str]

    planes = {}  # type: dict[str, cq.Face]
    partitions = {}  # type: dict[str, cq.Workplane]
    for face in faces:
        planes[face] = cq.Face.makePlane(
            basePnt=basePnt[face], dir=direction[face]
        )  # planes composing the Rve box
        partitions[face] = wk_plane.split(
            cq.Workplane().add(planes[face])
        )  # each partition cuts the object by the corresponding plane

        # detection of intersected faces
        if len(partitions[face].solids().all()) > 1:
            intersected_faces.append(face)

    if "x-" in intersected_faces and "x+" in intersected_faces:
        intersected_faces.remove("x-")
        intersected_faces.remove("x+")
        warnings.warn(
            "Object intersecting x+ and x- faces: not doing anything in this direction"
        )
    if "y-" in intersected_faces and "y+" in intersected_faces:
        intersected_faces.remove("y-")
        intersected_faces.remove("y+")
        warnings.warn(
            "Object intersecting y+ and y- faces: not doing anything in this direction"
        )
    if "z-" in intersected_faces and "z+" in intersected_faces:
        intersected_faces.remove("z-")
        intersected_faces.remove("z+")
        warnings.warn(
            "Object intersecting z+ and z- faces: not doing anything in this direction"
        )

    if len(intersected_faces) == 0:  # if no intersected faces = nothing to do
        periodic_object.append(wk_plane)

    elif len(intersected_faces) == 1:  # one face intersected
        f_0 = intersected_faces[0]
        periodic_object.append(
            partitions[f_0]
            .solids(face_dir[f_0])  # add the part of the
            .intersect(rve.box)  # object included in the rve
        )
        periodic_object.append(
            partitions[f_0]
            .solids(inverse_face_dir[f_0])  # translate the outside part of
            .translate(translate[f_0])  # the object in the rve and add
            .intersect(rve.box)  # it to the final object
        )

    elif len(intersected_faces) == 2:  # two faces intersected (edge)
        f_0 = intersected_faces[0]
        f_1 = intersected_faces[1]

        part = (
            cq.Workplane()
            .add(partitions[f_0].solids(face_dir[f_0]))
            .split(cq.Workplane().add(planes[f_1]))
        )
        periodic_object.append(part.solids(face_dir[f_1]).intersect(rve.box))
        periodic_object.append(
            part.solids(inverse_face_dir[f_1])
            .translate(translate[f_1])
            .intersect(rve.box)
        )

        part = (
            cq.Workplane()
            .add(partitions[f_0].solids(inverse_face_dir[f_0]))
            .split(cq.Workplane().add(planes[f_1]))
        )
        periodic_object.append(
            part.solids(face_dir[f_1]).translate(translate[f_0]).intersect(rve.box)
        )
        tslt = (
            translate[f_0][0] + translate[f_1][0],
            translate[f_0][1] + translate[f_1][1],
            translate[f_0][2] + translate[f_1][2],
        )
        periodic_object.append(
            part.solids(inverse_face_dir[f_1]).translate(tslt).intersect(rve.box)
        )

    elif len(intersected_faces) == 3:  # three faces intersected (corner)
        f_0 = intersected_faces[0]
        f_1 = intersected_faces[1]
        f_2 = intersected_faces[2]

        new_part = (
            cq.Workplane()
            .add(partitions[f_0].solids(face_dir[f_0]))
            .split(cq.Workplane().add(planes[f_1]))
            .solids(face_dir[f_1])
        )
        periodic_object.append(new_part.solids(face_dir[f_2]).intersect(rve.box))
        periodic_object.append(
            new_part.solids(inverse_face_dir[f_2])
            .translate(translate[f_2])
            .intersect(rve.box)
        )

        new_part = (
            cq.Workplane()
            .add(partitions[f_0].solids(face_dir[f_0]))
            .split(cq.Workplane().add(planes[f_1]))
            .solids(inverse_face_dir[f_1])
        )
        periodic_object.append(
            new_part.solids(face_dir[f_1]).translate(translate[f_1]).intersect(rve.box)
        )
        periodic_object.append(
            new_part.solids(inverse_face_dir[f_1])
            .translate((0, translate[f_1][1], translate[f_2][2]))
            .intersect(rve.box)
        )

        new_part = (
            cq.Workplane()
            .add(partitions[f_0].solids(inverse_face_dir[f_0]))
            .split(cq.Workplane().add(planes[f_1]))
            .solids(face_dir[f_1])
        )
        periodic_object.append(
            new_part.solids(face_dir[f_2]).translate(translate[f_0]).intersect(rve.box)
        )
        periodic_object.append(
            new_part.solids(inverse_face_dir[f_2])
            .translate((translate[f_0][0], 0, translate[f_2][2]))
            .intersect(rve.box)
        )

        new_part = (
            cq.Workplane()
            .add(partitions[f_0].solids(inverse_face_dir[f_0]))
            .split(cq.Workplane().add(planes[f_1]))
            .solids(inverse_face_dir[f_1])
        )
        periodic_object.append(
            new_part.solids(face_dir[f_2])
            .translate((translate[f_0][0], translate[f_1][1], 0))
            .intersect(rve.box)
        )
        periodic_object.append(
            new_part.solids(inverse_face_dir[f_2])
            .translate((translate[f_0][0], translate[f_1][1], translate[f_2][2]))
            .intersect(rve.box)
        )

    listSolids = [wp.val().Solids() for wp in periodic_object]
    flat_list = [solid.copy() for solids in listSolids for solid in solids]
    to_fuse = [cq.Shape(solid.wrapped) for solid in flat_list]
    shape = fuseShapes(cqShapeList=to_fuse, retain_edges=False)

    return Phase(shape=shape)
