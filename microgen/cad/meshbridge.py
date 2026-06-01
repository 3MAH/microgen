"""Mesh ↔ BREP bridges.

Functions in this module convert triangle meshes (numpy point/index arrays,
or pyvista ``PolyData``) into OCCT ``CadShape`` representations, or vice
versa (via the implicit ``Shape`` → BREP bridge in :func:`shape_to_cad`).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from ._install import require_cad
from .shape import CadShape, ShellCreationError, _topods_cast

if TYPE_CHECKING:
    from collections.abc import Sequence

    from ..shape._types import BoundsType
    from ..shape.shape import Shape


def mesh_to_shape(
    points: npt.NDArray[np.float64],
    triangles: npt.NDArray[np.int64],
) -> CadShape:
    """Convert a triangle mesh to a :class:`CadShape` via ``Poly_Triangulation``.

    One ``TopoDS_Face`` carries the full triangulation (OCCT native tessellated
    BREP representation).  This is the SOTA fast path: O(N) in pure OCCT C++
    with no Python-per-triangle overhead, and exports cleanly to STEP AP242
    (tessellated) and STL.

    :param points: ``(N, 3)`` array of vertex coordinates
    :param triangles: ``(M, 3)`` array of 0-indexed triangle vertex indices
    :return: wrapped ``TopoDS_Shell`` containing one tessellated face
    :raises ShellCreationError: if the triangulation cannot be built
    """
    require_cad()
    from OCP.BRep import BRep_Builder  # noqa: PLC0415
    from OCP.gp import gp_Pnt  # noqa: PLC0415
    from OCP.Poly import Poly_Triangle, Poly_Triangulation  # noqa: PLC0415
    from OCP.TopoDS import TopoDS_Face, TopoDS_Shell  # noqa: PLC0415

    pts = np.asarray(points, dtype=np.float64)
    tris = np.asarray(triangles, dtype=np.int64)
    if pts.ndim != 2 or pts.shape[1] != 3:
        err_msg = f"points must be (N, 3), got {pts.shape}"
        raise ValueError(err_msg)
    if tris.ndim != 2 or tris.shape[1] != 3:
        err_msg = f"triangles must be (M, 3), got {tris.shape}"
        raise ValueError(err_msg)
    if tris.size == 0:
        err_msg = "Cannot build a shell from an empty triangle list"
        raise ShellCreationError(err_msg)

    nb_nodes = int(pts.shape[0])
    nb_tri = int(tris.shape[0])
    triangulation = Poly_Triangulation(nb_nodes, nb_tri, False)
    for i in range(nb_nodes):
        triangulation.SetNode(
            i + 1, gp_Pnt(float(pts[i, 0]), float(pts[i, 1]), float(pts[i, 2]))
        )
    for i in range(nb_tri):
        a, b, c = int(tris[i, 0]), int(tris[i, 1]), int(tris[i, 2])
        triangulation.SetTriangle(i + 1, Poly_Triangle(a + 1, b + 1, c + 1))

    builder = BRep_Builder()
    face = TopoDS_Face()
    try:
        builder.MakeFace(face, triangulation)
    except Exception as err:
        err_msg = "OCCT refused the triangulation — check bounds and field."
        raise ShellCreationError(err_msg) from err

    shell = TopoDS_Shell()
    builder.MakeShell(shell)
    builder.Add(shell, face)
    return CadShape(shell)


def shape_to_cad(
    shape: Shape,
    bounds: BoundsType | None = None,
    resolution: int = 50,
) -> CadShape:
    """Bridge an implicit :class:`~microgen.shape.shape.Shape` to a CAD BREP.

    Runs the shape's :meth:`generate_surface_mesh` (marching cubes on the
    SDF for free Shapes; native renderer for concrete subclasses) then
    wraps the resulting triangle mesh into a single tessellated BREP face
    via :func:`mesh_to_shape`.

    Concrete subclasses with native primitive paths (``Box``, ``Sphere``,
    ``Tpms``, ``Spinodoid`` …) usually skip this bridge and call their
    own ``generate_cad``.  This function is the generic fallback for
    bare ``Shape`` instances built via :func:`microgen.shape.implicit_ops.from_field`
    or boolean composition (``a | b``, ``a - b``, ...).

    Requires the optional ``[cad]`` install extra.

    :param shape: the implicit shape to materialise
    :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)``; defaults to
        ``shape.bounds`` if set, else raises ``ValueError``
    :param resolution: marching-cubes grid resolution per axis
    :return: :class:`CadShape` wrapping the tessellated ``TopoDS_Shell``
    """
    require_cad()
    if shape.func is None:
        err_msg = "No implicit field defined — cannot build BREP from an empty Shape"
        raise NotImplementedError(err_msg)

    mesh = shape.generate_surface_mesh(bounds=bounds, resolution=resolution)
    if mesh.n_cells == 0:
        err_msg = "Generated mesh is empty — check bounds and field function"
        raise ValueError(err_msg)

    if not mesh.is_all_triangles:
        mesh.triangulate(inplace=True)
    triangles = mesh.faces.reshape(-1, 4)[:, 1:]
    points = np.asarray(mesh.points, dtype=np.float64)

    from ..shape.shape import ShellCreationError as _ShapeShellError  # noqa: PLC0415

    try:
        return mesh_to_shape(points, triangles)
    except Exception as err:
        err_msg = (
            "Failed to build the OCCT shell from the mesh; "
            "try to increase the resolution or adjust bounds."
        )
        raise _ShapeShellError(err_msg) from err


def _triangle_components(
    triangles: npt.NDArray[np.int64],
) -> npt.NDArray[np.int64]:
    """Label connected components of a triangle mesh by edge adjacency."""
    from collections import defaultdict  # noqa: PLC0415

    n_tri = int(triangles.shape[0])
    edges_sorted = np.sort(
        np.vstack(
            [triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]],
        ),
        axis=1,
    )
    _keys, inv = np.unique(edges_sorted, axis=0, return_inverse=True)
    owner = np.tile(np.arange(n_tri), 3)

    edge_to_tris: dict[int, list[int]] = defaultdict(list)
    for global_i, key_i in enumerate(inv):
        edge_to_tris[int(key_i)].append(int(owner[global_i]))

    adj: list[list[int]] = [[] for _ in range(n_tri)]
    for tlist in edge_to_tris.values():
        if len(tlist) == 2:
            adj[tlist[0]].append(tlist[1])
            adj[tlist[1]].append(tlist[0])

    component = -np.ones(n_tri, dtype=np.int64)
    n_comp = 0
    for start in range(n_tri):
        if component[start] >= 0:
            continue
        component[start] = n_comp
        stack = [start]
        while stack:
            t = stack.pop()
            for nb in adj[t]:
                if component[nb] < 0:
                    component[nb] = n_comp
                    stack.append(nb)
        n_comp += 1
    return component


def _walk_boundary_loops(
    comp_tris: npt.NDArray[np.int64],
) -> list[list[int]]:
    """Extract directed closed loops from the boundary edges of a triangle group.

    A boundary edge appears in exactly one triangle of the group; the loops
    are the closed chains of those edges in their original triangle direction.
    """
    from collections import defaultdict  # noqa: PLC0415

    ce = np.vstack(
        [comp_tris[:, [0, 1]], comp_tris[:, [1, 2]], comp_tris[:, [2, 0]]],
    )
    cs = np.sort(ce, axis=1)
    _keys, inv, counts = np.unique(
        cs,
        axis=0,
        return_inverse=True,
        return_counts=True,
    )
    bedges = ce[counts[inv] == 1]
    if len(bedges) == 0:
        return []

    used = np.zeros(len(bedges), dtype=bool)
    start_map: dict[int, list[int]] = defaultdict(list)
    for i, (a, _b) in enumerate(bedges):
        start_map[int(a)].append(i)

    loops: list[list[int]] = []
    for start_i in range(len(bedges)):
        if used[start_i]:
            continue
        cur = start_i
        used[cur] = True
        loop = [int(bedges[cur, 0])]
        while True:
            b = int(bedges[cur, 1])
            loop.append(b)
            if b == loop[0]:
                break
            nxt = next(
                (cand for cand in start_map[b] if not used[cand]),
                None,
            )
            if nxt is None:
                break
            used[nxt] = True
            cur = nxt
        if len(loop) >= 4 and loop[-1] == loop[0]:
            loops.append(loop)
    return loops


def mesh_to_planar_face(  # noqa: C901
    points: npt.NDArray[np.float64],
    triangles: npt.NDArray[np.int64],
    plane_origin: Sequence[float],
    plane_normal: Sequence[float],
) -> CadShape:
    """Build planar BREP face(s) whose wires trace the triangle-group boundary.

    Each connected component of the triangle group becomes a
    :class:`TopoDS_Face` whose underlying surface is ``Geom_Plane`` and
    whose outer/inner wires follow the boundary edges of the group on the
    plane.  Suitable both for STEP export (the BRep represents the right
    region, not a bounding rectangle) and for gmsh ``setPeriodic`` (the
    plane equation is recognised and the trimmed wires define the slave/
    master mesh region identically on opposite cell sides).

    All triangle vertices must lie on the plane within OCCT tolerance.

    :param points: ``(N, 3)`` array of vertex coordinates
    :param triangles: ``(M, 3)`` array of 0-indexed triangle vertex indices
    :param plane_origin: a point on the plane
    :param plane_normal: the plane's outward unit normal
    :return: wrapped face (single component) or shell (multiple components)
    :raises ShellCreationError: if no usable boundary loop is found
    """
    require_cad()
    from OCP.BRep import BRep_Builder  # noqa: PLC0415
    from OCP.BRepBuilderAPI import (  # noqa: PLC0415
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeWire,
    )
    from OCP.gp import gp_Dir, gp_Pln, gp_Pnt  # noqa: PLC0415
    from OCP.TopoDS import TopoDS_Shell  # noqa: PLC0415

    cast_wire = _topods_cast("Wire")

    pts = np.asarray(points, dtype=np.float64)
    tris = np.asarray(triangles, dtype=np.int64)
    if tris.size == 0:
        err_msg = "Cannot build a planar face from an empty triangle list"
        raise ShellCreationError(err_msg)

    plane = gp_Pln(
        gp_Pnt(float(plane_origin[0]), float(plane_origin[1]), float(plane_origin[2])),
        gp_Dir(float(plane_normal[0]), float(plane_normal[1]), float(plane_normal[2])),
    )

    n = np.asarray(plane_normal, dtype=np.float64)
    n = n / (float(np.linalg.norm(n)) or 1.0)
    helper = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u_ax = np.cross(n, helper)
    u_ax /= float(np.linalg.norm(u_ax)) or 1.0
    v_ax = np.cross(n, u_ax)
    o_arr = np.asarray(plane_origin, dtype=np.float64)

    component = _triangle_components(tris)
    n_comp = int(component.max()) + 1 if component.size else 0

    pnt_cache: dict[int, gp_Pnt] = {}

    def _occ_pnt(idx: int) -> gp_Pnt:
        cached = pnt_cache.get(idx)
        if cached is None:
            v = pts[idx]
            cached = gp_Pnt(float(v[0]), float(v[1]), float(v[2]))
            pnt_cache[idx] = cached
        return cached

    def _signed_area(loop: list[int]) -> float:
        rel = pts[loop[:-1]] - o_arr
        u = rel @ u_ax
        v = rel @ v_ax
        return 0.5 * float(np.sum(u * np.roll(v, -1) - np.roll(u, -1) * v))

    def _build_wire(loop: list[int]):
        wb = BRepBuilderAPI_MakeWire()
        for k in range(len(loop) - 1):
            edge_builder = BRepBuilderAPI_MakeEdge(
                _occ_pnt(loop[k]), _occ_pnt(loop[k + 1])
            )
            if not edge_builder.IsDone():
                err_msg = "BRepBuilderAPI_MakeEdge failed for boundary segment"
                raise ShellCreationError(err_msg)
            wb.Add(edge_builder.Edge())
        if not wb.IsDone():
            err_msg = "BRepBuilderAPI_MakeWire failed for boundary loop"
            raise ShellCreationError(err_msg)
        return wb.Wire()

    def _build_component_face(loops: list[list[int]]):
        areas = [_signed_area(L) for L in loops]
        outer_local = int(np.argmax([abs(a) for a in areas]))
        if areas[outer_local] < 0:
            loops = [list(reversed(L)) for L in loops]
        wires = [_build_wire(L) for L in loops]
        face_builder = BRepBuilderAPI_MakeFace(plane, wires[outer_local])
        for k, w in enumerate(wires):
            if k == outer_local:
                continue
            face_builder.Add(cast_wire(w.Reversed()))
        if not face_builder.IsDone():
            err_msg = "BRepBuilderAPI_MakeFace failed building a planar face"
            raise ShellCreationError(err_msg)
        return face_builder.Face()

    component_faces = []
    for ci in range(n_comp):
        comp_tris = tris[np.where(component == ci)[0]]
        loops = _walk_boundary_loops(comp_tris)
        if loops:
            component_faces.append(_build_component_face(loops))

    if not component_faces:
        err_msg = "No planar face could be built from the triangle group"
        raise ShellCreationError(err_msg)

    if len(component_faces) == 1:
        return CadShape(component_faces[0])

    builder = BRep_Builder()
    shell = TopoDS_Shell()
    builder.MakeShell(shell)
    for f in component_faces:
        builder.Add(shell, f)
    return CadShape(shell)


def mesh_to_shell_brep(
    points: npt.NDArray[np.float64],
    triangles: npt.NDArray[np.int64],
) -> CadShape:
    """Convert a triangle mesh to a shell with one planar BREP face per triangle.

    Slower than :func:`mesh_to_shape` (which uses a single tessellated face)
    but produces a shell with *real* geometric surfaces — required whenever
    the resulting shell is used as a cutting tool in boolean ops
    (``BRepAlgoAPI_Cut``, ``Workplane.split()``, etc.), which refuse a
    tessellated-only face.

    :param points: ``(N, 3)`` array of vertex coordinates
    :param triangles: ``(M, 3)`` array of 0-indexed triangle vertex indices
    :return: wrapped ``TopoDS_Shell`` with one planar face per triangle
    :raises ShellCreationError: if any triangle cannot be built into a face
    """
    require_cad()
    from OCP.BRep import BRep_Builder  # noqa: PLC0415
    from OCP.BRepBuilderAPI import (  # noqa: PLC0415
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeWire,
    )
    from OCP.gp import gp_Pnt  # noqa: PLC0415
    from OCP.TopoDS import TopoDS_Shell  # noqa: PLC0415

    pts = np.asarray(points, dtype=np.float64)
    tris = np.asarray(triangles, dtype=np.int64)
    if tris.size == 0:
        err_msg = "Cannot build a shell from an empty triangle list"
        raise ShellCreationError(err_msg)

    occ_points = [gp_Pnt(float(p[0]), float(p[1]), float(p[2])) for p in pts]

    builder = BRep_Builder()
    shell = TopoDS_Shell()
    builder.MakeShell(shell)

    try:
        for a, b, c in tris:
            e1 = BRepBuilderAPI_MakeEdge(occ_points[int(a)], occ_points[int(b)]).Edge()
            e2 = BRepBuilderAPI_MakeEdge(occ_points[int(b)], occ_points[int(c)]).Edge()
            e3 = BRepBuilderAPI_MakeEdge(occ_points[int(c)], occ_points[int(a)]).Edge()
            wire = BRepBuilderAPI_MakeWire(e1, e2, e3).Wire()
            face = BRepBuilderAPI_MakeFace(wire).Face()
            builder.Add(shell, face)
    except Exception as err:
        err_msg = (
            "Failed to build the OCCT shell from the mesh; "
            "try to increase the resolution or adjust bounds."
        )
        raise ShellCreationError(err_msg) from err

    return CadShape(shell)


def mesh_to_sewn_shell(
    points: npt.NDArray[np.float64],
    triangles: npt.NDArray[np.int64],
    tolerance: float | None = None,
) -> CadShape:
    """Convert a triangle mesh to a sewn shell with shared edges.

    Builds one planar BREP face per triangle then sews them via
    :class:`BRepBuilderAPI_Sewing` so coincident edges/vertices are merged
    into shared topology. The resulting shell has valid topology, which
    makes any subsequent boolean op (``Cut``, ``Common``, ``Fuse``) tractable
    — unlike :func:`mesh_to_shell_brep` whose triangle-soup output is
    pathologically slow for booleans.

    :param points: ``(N, 3)`` array of vertex coordinates
    :param triangles: ``(M, 3)`` array of 0-indexed triangle vertex indices
    :param tolerance: sewing tolerance; defaults to ``1e-6 * bbox_diag``
    :return: wrapped sewn ``TopoDS_Shape`` (Shell, or Compound of shells if
        the input is disconnected)
    :raises ShellCreationError: if any triangle cannot be built into a face
    """
    require_cad()
    from OCP.BRepBuilderAPI import (  # noqa: PLC0415
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeWire,
        BRepBuilderAPI_Sewing,
    )
    from OCP.gp import gp_Pnt  # noqa: PLC0415

    pts = np.asarray(points, dtype=np.float64)
    tris = np.asarray(triangles, dtype=np.int64)
    if tris.size == 0:
        err_msg = "Cannot build a shell from an empty triangle list"
        raise ShellCreationError(err_msg)

    if tolerance is None:
        bbox_diag = float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))
        tolerance = max(1e-9, 1e-6 * bbox_diag)

    occ_points = [gp_Pnt(float(p[0]), float(p[1]), float(p[2])) for p in pts]

    sewing = BRepBuilderAPI_Sewing(tolerance)
    try:
        for a, b, c in tris:
            e1 = BRepBuilderAPI_MakeEdge(occ_points[int(a)], occ_points[int(b)]).Edge()
            e2 = BRepBuilderAPI_MakeEdge(occ_points[int(b)], occ_points[int(c)]).Edge()
            e3 = BRepBuilderAPI_MakeEdge(occ_points[int(c)], occ_points[int(a)]).Edge()
            wire = BRepBuilderAPI_MakeWire(e1, e2, e3).Wire()
            face = BRepBuilderAPI_MakeFace(wire).Face()
            sewing.Add(face)
    except Exception as err:
        err_msg = (
            "Failed to build the OCCT shell from the mesh; "
            "try to increase the resolution or adjust bounds."
        )
        raise ShellCreationError(err_msg) from err

    sewing.Perform()
    return CadShape(sewing.SewedShape())
