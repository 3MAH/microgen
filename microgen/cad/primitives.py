"""OCCT primitive builders.

Free factory functions producing :class:`CadShape` instances for the core
parametric primitives (box, sphere, cylinder, capsule, ellipsoid,
polyhedron, extruded polygon).  These are the back-ends the implicit shape
classes (``Box``, ``Sphere``, …) call from their ``generate_cad`` methods.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ._install import require_cad
from .shape import CadShape, ShellCreationError, _topods_cast

if TYPE_CHECKING:
    from collections.abc import Sequence


def make_box(dim: Sequence[float], center: Sequence[float]) -> CadShape:
    """Axis-aligned box of size ``dim`` centered at ``center``."""
    require_cad()
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox  # noqa: PLC0415
    from OCP.gp import gp_Pnt  # noqa: PLC0415

    dx, dy, dz = (float(d) for d in dim)
    cx, cy, cz = (float(c) for c in center)
    corner = gp_Pnt(cx - dx / 2.0, cy - dy / 2.0, cz - dz / 2.0)
    return CadShape(BRepPrimAPI_MakeBox(corner, dx, dy, dz).Shape())


def make_sphere(radius: float, center: Sequence[float]) -> CadShape:
    """Sphere of given radius at ``center``."""
    require_cad()
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere  # noqa: PLC0415
    from OCP.gp import gp_Pnt  # noqa: PLC0415

    pnt = gp_Pnt(float(center[0]), float(center[1]), float(center[2]))
    return CadShape(BRepPrimAPI_MakeSphere(pnt, float(radius)).Shape())


def make_cylinder(
    radius: float,
    height: float,
    center: Sequence[float],
    axis: Sequence[float] = (1.0, 0.0, 0.0),
) -> CadShape:
    """Cylinder of given radius and height centered at ``center`` along ``axis``."""
    require_cad()
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder  # noqa: PLC0415
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt  # noqa: PLC0415

    h = float(height)
    ax_vec = np.asarray(axis, dtype=np.float64)
    ax_vec = ax_vec / np.linalg.norm(ax_vec)
    base = (
        float(center[0]) - h / 2.0 * ax_vec[0],
        float(center[1]) - h / 2.0 * ax_vec[1],
        float(center[2]) - h / 2.0 * ax_vec[2],
    )
    ax = gp_Ax2(
        gp_Pnt(*base),
        gp_Dir(float(ax_vec[0]), float(ax_vec[1]), float(ax_vec[2])),
    )
    return CadShape(BRepPrimAPI_MakeCylinder(ax, float(radius), h).Shape())


def make_capsule(
    radius: float,
    height: float,
    center: Sequence[float],
) -> CadShape:
    """Capsule (cylinder along X with hemispherical caps)."""
    require_cad()
    from OCP.BRepPrimAPI import (  # noqa: PLC0415
        BRepPrimAPI_MakeCylinder,
        BRepPrimAPI_MakeSphere,
    )
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt  # noqa: PLC0415

    cx, cy, cz = (float(c) for c in center)
    h = float(height)
    r = float(radius)
    base_axis = gp_Ax2(gp_Pnt(cx - h / 2.0, cy, cz), gp_Dir(1.0, 0.0, 0.0))
    cyl = BRepPrimAPI_MakeCylinder(base_axis, r, h).Shape()
    left = BRepPrimAPI_MakeSphere(gp_Pnt(cx - h / 2.0, cy, cz), r).Shape()
    right = BRepPrimAPI_MakeSphere(gp_Pnt(cx + h / 2.0, cy, cz), r).Shape()
    return CadShape(cyl).fuse(CadShape(left)).fuse(CadShape(right))


def make_ellipsoid(radii: Sequence[float], center: Sequence[float]) -> CadShape:
    """Ellipsoid of the given axis-aligned radii at ``center``.

    Built as a unit sphere transformed by a non-uniform scaling matrix.
    """
    require_cad()
    from OCP.BRepBuilderAPI import BRepBuilderAPI_GTransform  # noqa: PLC0415
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere  # noqa: PLC0415
    from OCP.gp import gp_GTrsf, gp_Mat, gp_Pnt, gp_XYZ  # noqa: PLC0415

    rx, ry, rz = (float(r) for r in radii)
    cx, cy, cz = (float(c) for c in center)
    sphere = BRepPrimAPI_MakeSphere(gp_Pnt(0.0, 0.0, 0.0), 1.0).Shape()

    gtrsf = gp_GTrsf()
    gtrsf.SetVectorialPart(
        gp_Mat(
            rx,
            0.0,
            0.0,
            0.0,
            ry,
            0.0,
            0.0,
            0.0,
            rz,
        ),
    )
    gtrsf.SetTranslationPart(gp_XYZ(cx, cy, cz))
    return CadShape(BRepBuilderAPI_GTransform(sphere, gtrsf, True).Shape())


def make_polyhedron(
    vertices: Sequence[Sequence[float]],
    faces_ixs: Sequence[Sequence[int]],
    center: Sequence[float] = (0.0, 0.0, 0.0),
) -> CadShape:
    """Polyhedron from vertex list and face→vertex index list.

    Each face's vertex index list must be closed (last == first).

    Uses ``BRepBuilderAPI_Sewing`` to stitch the independently-constructed
    faces into a shared-edge shell with consistent outward orientation,
    then ``BRepBuilderAPI_MakeSolid`` to close it.  (Raw ``BRep_Builder.Add``
    does not orient; the resulting solid would have mixed-sign volume.)
    """
    require_cad()
    from OCP.BRepBuilderAPI import (  # noqa: PLC0415
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeSolid,
        BRepBuilderAPI_MakeWire,
        BRepBuilderAPI_Sewing,
    )
    from OCP.gp import gp_Pnt  # noqa: PLC0415
    from OCP.ShapeFix import ShapeFix_Solid  # noqa: PLC0415
    from OCP.TopAbs import TopAbs_SHELL  # noqa: PLC0415
    from OCP.TopExp import TopExp_Explorer  # noqa: PLC0415

    cx, cy, cz = (float(c) for c in center)
    points = [
        gp_Pnt(float(v[0]) + cx, float(v[1]) + cy, float(v[2]) + cz) for v in vertices
    ]

    sewing = BRepBuilderAPI_Sewing()
    for ixs in faces_ixs:
        wire_builder = BRepBuilderAPI_MakeWire()
        for i1, i2 in zip(ixs, ixs[1:], strict=False):
            edge_builder = BRepBuilderAPI_MakeEdge(points[i1], points[i2])
            if not edge_builder.IsDone():
                err_msg = "BRepBuilderAPI_MakeEdge failed for polyhedron face"
                raise ShellCreationError(err_msg)
            wire_builder.Add(edge_builder.Edge())
        if not wire_builder.IsDone():
            err_msg = "BRepBuilderAPI_MakeWire failed for polyhedron face"
            raise ShellCreationError(err_msg)
        face_builder = BRepBuilderAPI_MakeFace(wire_builder.Wire())
        if not face_builder.IsDone():
            err_msg = "BRepBuilderAPI_MakeFace failed for polyhedron face"
            raise ShellCreationError(err_msg)
        sewing.Add(face_builder.Face())
    sewing.Perform()
    sewn = sewing.SewedShape()

    # Extract the shell (sewing may return it directly or wrapped in a compound).
    exp = TopExp_Explorer(sewn, TopAbs_SHELL)
    if not exp.More():
        err_msg = "Sewing did not produce a shell — check face connectivity"
        raise ShellCreationError(err_msg)
    shell = _topods_cast("Shell")(exp.Current())

    solid_builder = BRepBuilderAPI_MakeSolid(shell)
    if not solid_builder.IsDone():
        err_msg = "BRepBuilderAPI_MakeSolid failed; sewn shell is not closed"
        raise ShellCreationError(err_msg)
    fixer = ShapeFix_Solid(solid_builder.Solid())
    fixer.Perform()
    return CadShape(fixer.Solid())


def make_extruded_polygon(
    list_corners: Sequence[tuple[float, float]],
    height: float,
    center: Sequence[float],
) -> CadShape:
    """Extrude a 2D polygon (in the YZ plane) along the X axis."""
    require_cad()
    from OCP.BRepBuilderAPI import (  # noqa: PLC0415
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeWire,
    )
    from OCP.BRepPrimAPI import BRepPrimAPI_MakePrism  # noqa: PLC0415
    from OCP.gp import gp_Pnt, gp_Vec  # noqa: PLC0415

    cx, cy, cz = (float(c) for c in center)
    h = float(height)
    x_base = cx - h / 2.0
    pts = [gp_Pnt(x_base, cy + float(y), cz + float(z)) for (y, z) in list_corners]
    if (list_corners[0][0], list_corners[0][1]) != (
        list_corners[-1][0],
        list_corners[-1][1],
    ):
        pts.append(pts[0])

    wire_builder = BRepBuilderAPI_MakeWire()
    for i in range(len(pts) - 1):
        edge_builder = BRepBuilderAPI_MakeEdge(pts[i], pts[i + 1])
        if not edge_builder.IsDone():
            err_msg = "BRepBuilderAPI_MakeEdge failed for extruded polygon"
            raise ShellCreationError(err_msg)
        wire_builder.Add(edge_builder.Edge())
    if not wire_builder.IsDone():
        err_msg = "BRepBuilderAPI_MakeWire failed for extruded polygon"
        raise ShellCreationError(err_msg)
    face_builder = BRepBuilderAPI_MakeFace(wire_builder.Wire())
    if not face_builder.IsDone():
        err_msg = "BRepBuilderAPI_MakeFace failed for extruded polygon"
        raise ShellCreationError(err_msg)
    extruded = BRepPrimAPI_MakePrism(face_builder.Face(), gp_Vec(h, 0.0, 0.0)).Shape()
    return CadShape(extruded)
