"""Periodic cut — rearrange a shape so it wraps around an RVE.

Pure OCP implementation (``cadquery-ocp-novtk``); no ``cadquery`` dependency.
Requires the ``[cad]`` extra.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

from .cad import (
    CadShape,
    enumerate_solids,
    intersect_solids_with_box,
    make_compound_from_solids,
    make_plane_face,
    select_solids_on_side,
    split_shape,
    translate_solid,
)
from .operations import fuse_shapes
from .phase import Phase

if TYPE_CHECKING:
    from .rve import Rve


# -- Face labels and directions ------------------------------------------------
_NO_INTERSECTION = 0
_FACE = 1
_EDGE = 2
_CORNER = 3

_FACES = ["x-", "x+", "y-", "y+", "z-", "z+"]

# Outward normal at each RVE face; used both to pick the base point of the
# cutting plane and to pick which side a partitioned solid belongs to.
_DIRECTION: dict[str, tuple[int, int, int]] = {
    "x-": (-1, 0, 0),
    "x+": (1, 0, 0),
    "y-": (0, -1, 0),
    "y+": (0, 1, 0),
    "z-": (0, 0, -1),
    "z+": (0, 0, 1),
}

# For face "x-" (minus-X face of the RVE), ">X" = solids on the +X side of
# the cutting plane — i.e. inside the RVE.  "<X" = solids on the -X side —
# outside the RVE, to be translated back in.
_SIDE_DIR: dict[str, tuple[int, int, int]] = {
    "x-": (1, 0, 0),
    "x+": (-1, 0, 0),
    "y-": (0, 1, 0),
    "y+": (0, -1, 0),
    "z-": (0, 0, 1),
    "z+": (0, 0, -1),
}


def _base_pnt(face: str, rve: Rve) -> tuple[float, float, float]:
    if face == "x-":
        return (float(rve.min_point[0]), 0.0, 0.0)
    if face == "x+":
        return (float(rve.max_point[0]), 0.0, 0.0)
    if face == "y-":
        return (0.0, float(rve.min_point[1]), 0.0)
    if face == "y+":
        return (0.0, float(rve.max_point[1]), 0.0)
    if face == "z-":
        return (0.0, 0.0, float(rve.min_point[2]))
    if face == "z+":
        return (0.0, 0.0, float(rve.max_point[2]))
    err_msg = f"Unknown face label {face!r}"
    raise ValueError(err_msg)


def _translate(face: str, rve: Rve) -> tuple[float, float, float]:
    """Compute the wrap-around translation for solids on a face.

    The translation is applied to solids sitting on the 'outside' half of
    *face* so they re-enter the RVE through the opposite face.
    """
    if face == "x-":
        return (float(rve.dim[0]), 0.0, 0.0)
    if face == "x+":
        return (-float(rve.dim[0]), 0.0, 0.0)
    if face == "y-":
        return (0.0, float(rve.dim[1]), 0.0)
    if face == "y+":
        return (0.0, -float(rve.dim[1]), 0.0)
    if face == "z-":
        return (0.0, 0.0, float(rve.dim[2]))
    if face == "z+":
        return (0.0, 0.0, -float(rve.dim[2]))
    err_msg = f"Unknown face label {face!r}"
    raise ValueError(err_msg)


# -- Partitioning --------------------------------------------------------------


def _detect_intersected_faces(
    shape: CadShape,
    rve: Rve,
) -> tuple[list[str], dict[str, CadShape], dict[str, CadShape]]:
    """Split *shape* by each of the six RVE faces; report which faces cut it."""
    intersected_faces: list[str] = []
    rve_planes: dict[str, CadShape] = {}
    partitions: dict[str, CadShape] = {}

    for face in _FACES:
        plane = make_plane_face(_base_pnt(face, rve), _DIRECTION[face])
        rve_planes[face] = plane
        partitions[face] = split_shape(shape, plane)

        # A face is "intersected" iff the split produced more than one solid.
        if len(enumerate_solids(partitions[face])) > 1:
            intersected_faces.append(face)

    # Skip pairs of opposite faces (object straddles entire RVE in that axis).
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


# -- Intersection helpers ------------------------------------------------------


def _intersect_face(
    face: str,
    rve: Rve,
    partitions: dict[str, CadShape],
) -> list[CadShape]:
    """One-face intersection: half stays, half translates back via periodicity."""
    base = _base_pnt(face, rve)
    inside_solids = select_solids_on_side(partitions[face], base, _SIDE_DIR[face])
    outside_solids = select_solids_on_side(
        partitions[face],
        base,
        tuple(-d for d in _SIDE_DIR[face]),
    )
    tr = _translate(face, rve)
    translated = [translate_solid(s, tr) for s in outside_solids]

    return [
        intersect_solids_with_box(inside_solids, rve.box),
        intersect_solids_with_box(translated, rve.box),
    ]


def _intersect_edge(
    faces: tuple[str, str],
    rve: Rve,
    partitions: dict[str, CadShape],
    rve_planes: dict[str, CadShape],
) -> list[CadShape]:
    """Two-face (edge) intersection: 4 quadrants, 3 need translation."""
    f_0, f_1 = faces
    base_0 = _base_pnt(f_0, rve)
    base_1 = _base_pnt(f_1, rve)
    tr_0 = _translate(f_0, rve)
    tr_1 = _translate(f_1, rve)

    # (++) — inside both half-spaces, no translation
    p0_inside = select_solids_on_side(partitions[f_0], base_0, _SIDE_DIR[f_0])
    split_pp = split_shape(
        make_compound_from_solids(p0_inside),
        rve_planes[f_1],
    )
    pp = select_solids_on_side(split_pp, base_1, _SIDE_DIR[f_1])

    # (+-) -- inside first, outside second: translate by tr_1
    pm = select_solids_on_side(
        split_pp,
        base_1,
        tuple(-d for d in _SIDE_DIR[f_1]),
    )
    pm_t = [translate_solid(s, tr_1) for s in pm]

    # Split the "outside first" half by the second plane to get (-+) and (--)
    p0_outside = select_solids_on_side(
        partitions[f_0],
        base_0,
        tuple(-d for d in _SIDE_DIR[f_0]),
    )
    split_m = split_shape(
        make_compound_from_solids(p0_outside),
        rve_planes[f_1],
    )
    # (-+) — outside first, inside second: translate by tr_0
    mp = select_solids_on_side(split_m, base_1, _SIDE_DIR[f_1])
    mp_t = [translate_solid(s, tr_0) for s in mp]
    # (--) — outside both: translate by tr_0 + tr_1
    mm = select_solids_on_side(
        split_m,
        base_1,
        tuple(-d for d in _SIDE_DIR[f_1]),
    )
    tslt = (tr_0[0] + tr_1[0], tr_0[1] + tr_1[1], tr_0[2] + tr_1[2])
    mm_t = [translate_solid(s, tslt) for s in mm]

    return [
        intersect_solids_with_box(pp, rve.box),
        intersect_solids_with_box(pm_t, rve.box),
        intersect_solids_with_box(mp_t, rve.box),
        intersect_solids_with_box(mm_t, rve.box),
    ]


def _intersect_corner(
    faces: tuple[str, str, str],
    rve: Rve,
    partitions: dict[str, CadShape],
    rve_planes: dict[str, CadShape],
) -> list[CadShape]:
    """Three-face (corner) intersection: 8 octants, 7 need translation."""
    f_0, f_1, f_2 = faces
    base_0 = _base_pnt(f_0, rve)
    base_1 = _base_pnt(f_1, rve)
    base_2 = _base_pnt(f_2, rve)
    tr_0 = _translate(f_0, rve)
    tr_1 = _translate(f_1, rve)
    tr_2 = _translate(f_2, rve)

    def add(*vs: tuple[float, float, float]) -> tuple[float, float, float]:
        return (
            sum(v[0] for v in vs),
            sum(v[1] for v in vs),
            sum(v[2] for v in vs),
        )

    def neg(d: tuple[int, int, int]) -> tuple[int, int, int]:
        return (-d[0], -d[1], -d[2])

    results: list[CadShape] = []

    # Branch by sign of f_0 side
    for sign_0, tr_x in (
        (_SIDE_DIR[f_0], (0.0, 0.0, 0.0)),
        (neg(_SIDE_DIR[f_0]), tr_0),
    ):
        p0 = select_solids_on_side(partitions[f_0], base_0, sign_0)
        split_1 = split_shape(make_compound_from_solids(p0), rve_planes[f_1])

        # Branch by sign of f_1 side
        for sign_1, tr_y in (
            (_SIDE_DIR[f_1], (0.0, 0.0, 0.0)),
            (neg(_SIDE_DIR[f_1]), tr_1),
        ):
            p1 = select_solids_on_side(split_1, base_1, sign_1)
            split_2 = split_shape(make_compound_from_solids(p1), rve_planes[f_2])

            # Branch by sign of f_2 side
            for sign_2, tr_z in (
                (_SIDE_DIR[f_2], (0.0, 0.0, 0.0)),
                (neg(_SIDE_DIR[f_2]), tr_2),
            ):
                p2 = select_solids_on_side(split_2, base_2, sign_2)
                shift = add(tr_x, tr_y, tr_z)
                if shift == (0.0, 0.0, 0.0):
                    translated = p2
                else:
                    translated = [translate_solid(s, shift) for s in p2]
                results.append(intersect_solids_with_box(translated, rve.box))

    return results


# -- Public API ----------------------------------------------------------------


def periodic_split_and_translate(phase: Phase, rve: Rve) -> Phase:
    """Rearrange a phase periodically so it wraps the RVE.

    Pure OCP implementation — no cadquery dependency.

    :param phase: Phase to cut periodically
    :param rve: RVE for periodicity

    :return: resulting Phase
    """
    shape = phase.shape
    if shape is None:
        err_msg = "Cannot apply periodic_split_and_translate to an empty phase"
        raise ValueError(err_msg)

    intersected_faces, rve_planes, partitions = _detect_intersected_faces(shape, rve)

    periodic_object: list[CadShape] = []
    if len(intersected_faces) == _NO_INTERSECTION:
        periodic_object.append(shape)
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

    # Flatten: each CadShape in periodic_object wraps a compound of solids.
    # Fuse all resulting solids into one final CadShape.
    all_solids: list[Any] = []
    for piece in periodic_object:
        all_solids.extend(enumerate_solids(piece))

    if not all_solids:
        warnings.warn(
            "periodic_split_and_translate produced no solids",
            stacklevel=2,
        )
        return Phase(shape=shape)

    to_fuse = [CadShape(s) for s in all_solids]
    fused = fuse_shapes(to_fuse, retain_edges=False)
    return Phase(shape=fused)
