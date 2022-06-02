"""
Periodic function to cut a shape periodically according to a RVE
"""
import cadquery as cq

from .operations import fuseParts
from .rve import Rve


def periodic(cqshape: cq.Shape, rve: Rve) -> tuple:
    """
    Rearrange cqshape periodically according to the rve

    :param cqshape: CQ Shape to cut periodically
    :param rve: RVE for periodicity

    :return return_object_periodic[0].copy(): cutted object
    :return flat_list: list of cutted parts
    """

    wk_plane = cq.Workplane().add(cqshape.Solids())  # shape to cut
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
        "x-": (0, 0, 0),
        "x+": (rve.dx, 0, 0),
        "y-": (0, 0, 0),
        "y+": (0, rve.dy, 0),
        "z-": (0, 0, 0),
        "z+": (0, 0, rve.dz),
    }

    face_dir = {"x-": ">X", "x+": "<X", "y-": ">Y", "y+": "<Y", "z-": ">Z", "z+": "<Z"}
    inverse_face_dir = {"x-": "<X", "x+": ">X", "y-": "<Y", "y+": ">Y", "z-": "<Z", "z+": ">Z"}

    translate = {
        "x-": (rve.dx, 0, 0),
        "x+": (-rve.dx, 0, 0),
        "y-": (0, rve.dy, 0),
        "y+": (0, -rve.dy, 0),
        "z-": (0, 0, rve.dz),
        "z+": (0, 0, -rve.dz),
    }

    intersected_faces = []

    planes = {}
    partitions = {}
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

    # A REFLECHIR
    # Erreur/Warning quand un objet dépasse des deux côtés ?
    if "x-" in intersected_faces and "x+" in intersected_faces:
        intersected_faces.remove("x-")
        intersected_faces.remove("x+")
    if "y-" in intersected_faces and "y+" in intersected_faces:
        intersected_faces.remove("y-")
        intersected_faces.remove("y+")
    if "z-" in intersected_faces and "z+" in intersected_faces:
        intersected_faces.remove("z-")
        intersected_faces.remove("z+")

    if len(intersected_faces) == 0:  # if no intersected faces = nothing to do
        periodic_object.append(wk_plane)

    elif len(intersected_faces) == 1:  # one face intersected
        f_0 = intersected_faces[0]
        periodic_object.append(
            partitions[f_0]
            .solids(face_dir[f_0])  # add the part of the
            .intersect(rve.Box)  # object included in the rve
        )
        periodic_object.append(
            partitions[f_0]
            .solids(inverse_face_dir[f_0])  # translate the outside part of
            .translate(translate[f_0])  # the object in the rve and add
            .intersect(rve.Box)  # it to the final object
        )

    elif len(intersected_faces) == 2:  # two faces intersected (edge)
        f_0 = intersected_faces[0]
        f_1 = intersected_faces[1]

        part = (
            cq.Workplane()
            .add(partitions[f_0].solids(face_dir[f_0]))
            .split(cq.Workplane().add(planes[f_1]))
        )
        periodic_object.append(part.solids(face_dir[f_1]).intersect(rve.Box))
        periodic_object.append(
            part.solids(inverse_face_dir[f_1])
            .translate(translate[f_1])
            .intersect(rve.Box)
        )

        part = (
            cq.Workplane()
            .add(partitions[f_0].solids(inverse_face_dir[f_0]))
            .split(cq.Workplane().add(planes[f_1]))
        )
        periodic_object.append(
            part.solids(face_dir[f_1]).translate(translate[f_0]).intersect(rve.Box)
        )
        tslt = (
            translate[f_0][0] + translate[f_1][0],
            translate[f_0][1] + translate[f_1][1],
            translate[f_0][2] + translate[f_1][2],
        )
        periodic_object.append(
            part.solids(inverse_face_dir[f_1]).translate(tslt).intersect(rve.Box)
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
        periodic_object.append(new_part.solids(face_dir[f_2]).intersect(rve.Box))
        periodic_object.append(
            new_part.solids(inverse_face_dir[f_2])
            .translate(translate[f_2])
            .intersect(rve.Box)
        )

        new_part = (
            cq.Workplane()
            .add(partitions[f_0].solids(face_dir[f_0]))
            .split(cq.Workplane().add(planes[f_1]))
            .solids(inverse_face_dir[f_1])
        )
        periodic_object.append(
            new_part.solids(face_dir[f_1]).translate(translate[f_1]).intersect(rve.Box)
        )
        periodic_object.append(
            new_part.solids(inverse_face_dir[f_1])
            .translate((0, translate[f_1][1], translate[f_2][2]))
            .intersect(rve.Box)
        )

        new_part = (
            cq.Workplane()
            .add(partitions[f_0].solids(inverse_face_dir[f_0]))
            .split(cq.Workplane().add(planes[f_1]))
            .solids(face_dir[f_1])
        )
        periodic_object.append(
            new_part.solids(face_dir[f_2]).translate(translate[f_0]).intersect(rve.Box)
        )
        periodic_object.append(
            new_part.solids(inverse_face_dir[f_2])
            .translate((translate[f_0][0], 0, translate[f_2][2]))
            .intersect(rve.Box)
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
            .intersect(rve.Box)
        )
        periodic_object.append(
            new_part.solids(inverse_face_dir[f_2])
            .translate((translate[f_0][0], translate[f_1][1], translate[f_2][2]))
            .intersect(rve.Box)
        )

    occ_solids_list = [s.val().Solids() for s in periodic_object]
    flat_list = [item.copy() for sublist in occ_solids_list for item in sublist]
    to_fuse = [cq.Shape(s.wrapped) for s in flat_list]
    return_object_periodic = fuseParts(to_fuse, False)
    return (return_object_periodic[0].copy(), flat_list)
