"""Sew a triangulated periodic surface into a closed OCCT shell.

Given a marching-cubes (or otherwise periodic-aligned) triangle mesh of an
implicit field whose iso-surface meets the cell boundary on cap planes,
:func:`mesh_to_periodic_shell` groups boundary triangles per cap plane,
builds **one planar BREP face per cap** (via :func:`microgen.cad.mesh_to_planar_face`)
carrying the actual cap-wire trace, and sews them together with the
interior triangles into a closed shell suitable for boolean ops, STEP
export, and gmsh ``setPeriodic`` constraints.

Shared by :class:`microgen.shape.tpms.Tpms` and
:class:`microgen.shape.spinodoid.Spinodoid` — the two TPMS/F-rep shapes
that produce periodic meshes that need this cap-aware sewing.

Requires the optional ``[cad]`` install extra.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from collections.abc import Sequence

    from microgen.cad import CadShape

    from ._types import BoundsType


def mesh_to_periodic_shell(
    points: npt.NDArray[np.float64],
    triangles: npt.NDArray[np.int64],
    bounds: BoundsType | Sequence[float],
) -> CadShape:
    """Build a sewn OCCT shell from a triangulated surface in an AABB cell.

    Triangles whose three vertices share the mesh's exact extremum on an
    AABB cap plane are grouped per face and converted to a single planar
    BREP face via :func:`microgen.cad.mesh_to_planar_face` (so STEP shows
    the actual gyroid/spinodoid cuts, not a bounding cube, and gmsh
    ``setPeriodic`` can pair opposite cap faces by plane equation).
    Remaining triangles become one planar BREP face per triangle (raw,
    not pre-sewn — sewing must stitch them to the cap wires along the
    seam).  Everything is sewn into a closed shell via
    :class:`OCP.BRepBuilderAPI.BRepBuilderAPI_Sewing`.

    :param points: ``(N, 3)`` vertex coordinates
    :param triangles: ``(M, 3)`` vertex-index triplets
    :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)`` of the cell
    :return: :class:`microgen.cad.CadShape` wrapping the sewn ``TopoDS_Shell``
    """
    from OCP.BRepBuilderAPI import (  # noqa: PLC0415
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeWire,
        BRepBuilderAPI_Sewing,
    )
    from OCP.gp import gp_Pnt  # noqa: PLC0415
    from OCP.TopAbs import TopAbs_FACE  # noqa: PLC0415
    from OCP.TopExp import TopExp_Explorer  # noqa: PLC0415

    from microgen.cad import CadShape as _CadShape  # noqa: PLC0415
    from microgen.cad import mesh_to_planar_face  # noqa: PLC0415

    pts = np.asarray(points, dtype=np.float64)
    tris = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)

    drift_tol = 1e-9 * float(max(abs(b) for b in bounds) or 1.0)

    consumed = np.zeros(tris.shape[0], dtype=bool)
    on_plane: list[tuple[int, int, float, npt.NDArray[np.int64]]] = []
    for axis in range(3):
        for sign in (-1, +1):
            extremum = (
                float(pts[:, axis].max()) if sign > 0 else float(pts[:, axis].min())
            )
            expected = bounds[2 * axis + (1 if sign > 0 else 0)]
            if abs(extremum - expected) > drift_tol:
                continue
            vert_on = pts[:, axis] == extremum
            tri_on = np.all(vert_on[tris], axis=1) & ~consumed
            if tri_on.any():
                on_plane.append((axis, sign, extremum, np.where(tri_on)[0]))
                consumed |= tri_on
    interior_idx = np.where(~consumed)[0]

    bbox_diag = float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))
    sew_tol = max(1e-9, 1e-6 * bbox_diag)
    sewing = BRepBuilderAPI_Sewing(sew_tol)

    for axis, sign, extremum, tri_idx in on_plane:
        origin = [0.0, 0.0, 0.0]
        origin[axis] = extremum
        normal = [0.0, 0.0, 0.0]
        normal[axis] = float(sign)
        planar = mesh_to_planar_face(pts, tris[tri_idx], origin, normal)
        exp = TopExp_Explorer(planar.wrapped, TopAbs_FACE)
        while exp.More():
            sewing.Add(exp.Current())
            exp.Next()

    interior_tris = tris[interior_idx]
    used_vertices = np.unique(interior_tris)
    pnt_cache = {
        int(i): gp_Pnt(float(pts[i, 0]), float(pts[i, 1]), float(pts[i, 2]))
        for i in used_vertices
    }
    for a, b, c in interior_tris:
        ga, gb, gc = pnt_cache[int(a)], pnt_cache[int(b)], pnt_cache[int(c)]
        e1 = BRepBuilderAPI_MakeEdge(ga, gb).Edge()
        e2 = BRepBuilderAPI_MakeEdge(gb, gc).Edge()
        e3 = BRepBuilderAPI_MakeEdge(gc, ga).Edge()
        wire = BRepBuilderAPI_MakeWire(e1, e2, e3).Wire()
        face = BRepBuilderAPI_MakeFace(wire).Face()
        sewing.Add(face)

    sewing.Perform()
    return _CadShape(sewing.SewedShape())


__all__ = ["mesh_to_periodic_shell"]
