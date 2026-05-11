"""CAD backend — direct OCCT (via OCP) replacement for CadQuery.

=========================================================
CAD backend (:mod:`microgen.cad`)
=========================================================

All CadQuery calls in microgen have been replaced by direct OCP
(``cadquery-ocp-novtk``) calls housed in this module.  OCP is an *optional*
dependency — install via ``pip install microgen[cad]``.

The module's top-level body does not import OCP, so ``import microgen.cad``
always succeeds.  The OCP-dependent functions import it lazily and raise a
helpful ``ImportError`` with install instructions if OCP is missing.

Return type
-----------

CAD-producing functions return a :class:`CadShape` — a thin wrapper around
an OCCT ``TopoDS_Shape`` exposing ``.wrapped`` for downstream OCP calls, plus
convenience methods (``translate``, ``rotate``, ``fuse``, ``cut``,
``export_stl``, ``export_step``, ``export_brep``).  ``.wrapped`` matches the
attribute name CadQuery's ``Shape`` exposed, so most legacy call sites keep
working unchanged.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from pathlib import Path

    from OCP.TopoDS import TopoDS_Shape


_INSTALL_HINT = (
    "microgen's CAD backend requires the OCP (OCCT) Python bindings. "
    "Install with:  pip install 'microgen[cad]'  "
    "(this pulls cadquery-ocp-novtk; on conda-forge use `ocp` instead)."
)


class _Centre(tuple):
    """Tuple-like 3D point exposing ``.x``, ``.y``, ``.z`` and ``.toTuple()``.

    Returned by :meth:`CadShape.Center`.  Mimics just enough of
    ``cadquery.Vector`` to drop-in where existing code uses
    ``shape.Center().toTuple()`` or ``shape.Center().x``.
    """

    __slots__ = ()

    def __new__(cls, x: float, y: float, z: float) -> _Centre:
        """Create a 3-tuple ``(x, y, z)``."""
        return super().__new__(cls, (float(x), float(y), float(z)))

    @property
    def x(self) -> float:
        """X coordinate."""
        return self[0]

    @property
    def y(self) -> float:
        """Y coordinate."""
        return self[1]

    @property
    def z(self) -> float:
        """Z coordinate."""
        return self[2]

    def toTuple(self) -> tuple[float, float, float]:  # noqa: N802
        """CadQuery-compatible alias returning ``(x, y, z)``."""
        return (self[0], self[1], self[2])


class _BBox:
    """Axis-aligned bounding box exposing CadQuery-style ``xmin``/``xmax``/…

    Returned by :meth:`CadShape.BoundingBox`.  Also indexable as a 6-tuple
    ``(xmin, ymin, zmin, xmax, ymax, zmax)`` matching OCCT's ``Bnd_Box.Get``.
    """

    __slots__ = ("xmax", "xmin", "ymax", "ymin", "zmax", "zmin")

    def __init__(
        self,
        xmin: float,
        ymin: float,
        zmin: float,
        xmax: float,
        ymax: float,
        zmax: float,
    ) -> None:
        """Initialize from the 6 axis-aligned extents."""
        self.xmin = float(xmin)
        self.ymin = float(ymin)
        self.zmin = float(zmin)
        self.xmax = float(xmax)
        self.ymax = float(ymax)
        self.zmax = float(zmax)

    @property
    def DiagonalLength(self) -> float:  # noqa: N802
        """CadQuery-compatible diagonal length."""
        dx = self.xmax - self.xmin
        dy = self.ymax - self.ymin
        dz = self.zmax - self.zmin
        return float((dx * dx + dy * dy + dz * dz) ** 0.5)


def require_cad() -> None:
    """Raise :class:`ImportError` if the CAD backend (OCP) is not importable."""
    try:
        import OCP  # noqa: F401
    except ImportError as err:
        raise ImportError(_INSTALL_HINT) from err


def _run_boolean(op_cls: Any, a: CadShape, b: CadShape, label: str) -> TopoDS_Shape:
    """Run an OCCT boolean op and raise on failure.

    Older OCP releases don't expose ``HasErrors()`` / ``IsDone()`` on the
    ``BRepAlgoAPI_*`` classes; we probe via ``getattr`` and skip the check
    when the API isn't available.
    """
    op = op_cls(a.wrapped, b.wrapped)
    has_errors = getattr(op, "HasErrors", None)
    if callable(has_errors) and has_errors():
        err_msg = f"BRepAlgoAPI_{label} failed"
        raise RuntimeError(err_msg)
    return op.Shape()


def _topods_cast(name: str) -> Any:
    """Return ``TopoDS.<name>`` cast helper, tolerant of OCP version drift.

    Older OCP releases expose the static cast as ``TopoDS.Shell_s`` (pybind11
    ``_s`` convention); newer releases expose the unsuffixed ``TopoDS.Shell``.
    Try the suffixed form first, fall back to unsuffixed.
    """
    from OCP.TopoDS import TopoDS

    return getattr(TopoDS, f"{name}_s", None) or getattr(TopoDS, name)


# ---------------------------------------------------------------------------
# CadShape wrapper
# ---------------------------------------------------------------------------


class CadShape:
    """Thin wrapper around an OCCT ``TopoDS_Shape``.

    Preserves the ``.wrapped`` attribute name used by CadQuery so downstream
    OCP calls (``BRepAlgoAPI_Fuse(a.wrapped, b.wrapped)``) keep working.

    ``_mesh_volume`` (optional) is a trusted volume in the source mesh's
    units, set by mesh-derived constructors (e.g. the TPMS periodic shell)
    where OCCT's surface-integral volume is unreliable on invalid topology.
    :meth:`Volume` prefers it over the OCCT integral when present.
    """

    __slots__ = ("_mesh_volume", "wrapped")

    def __init__(self, shape: TopoDS_Shape) -> None:
        """Wrap an OCCT ``TopoDS_Shape``."""
        self.wrapped = shape
        self._mesh_volume: float | None = None

    # -- transforms --------------------------------------------------------

    def translate(self, offset: Sequence[float]) -> CadShape:
        """Return a translated copy."""
        from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
        from OCP.gp import gp_Trsf, gp_Vec

        trsf = gp_Trsf()
        trsf.SetTranslation(
            gp_Vec(float(offset[0]), float(offset[1]), float(offset[2]))
        )
        transformed = BRepBuilderAPI_Transform(self.wrapped, trsf, True).Shape()
        return CadShape(transformed)

    def rotate(
        self,
        center: Sequence[float],
        axis: Sequence[float],
        angle_degrees: float,
    ) -> CadShape:
        """Return a rotated copy (angle in degrees, axis is a unit vector)."""
        from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
        from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf

        trsf = gp_Trsf()
        ax = gp_Ax1(
            gp_Pnt(float(center[0]), float(center[1]), float(center[2])),
            gp_Dir(float(axis[0]), float(axis[1]), float(axis[2])),
        )
        trsf.SetRotation(ax, float(np.deg2rad(angle_degrees)))
        transformed = BRepBuilderAPI_Transform(self.wrapped, trsf, True).Shape()
        return CadShape(transformed)

    def copy(self) -> CadShape:
        """Return an independent copy (deep topology copy)."""
        from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy

        return CadShape(BRepBuilderAPI_Copy(self.wrapped).Shape())

    # -- boolean ops -------------------------------------------------------

    def fuse(self, other: CadShape) -> CadShape:
        """Boolean fusion: ``self ∪ other``."""
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse

        return CadShape(_run_boolean(BRepAlgoAPI_Fuse, self, other, "Fuse"))

    def cut(self, other: CadShape) -> CadShape:
        """Boolean difference: ``self \\ other``."""
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut

        return CadShape(_run_boolean(BRepAlgoAPI_Cut, self, other, "Cut"))

    def intersect(self, other: CadShape) -> CadShape:
        """Boolean intersection: ``self ∩ other``."""
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Common

        return CadShape(_run_boolean(BRepAlgoAPI_Common, self, other, "Common"))

    # -- topology queries --------------------------------------------------

    def solids(self) -> list[CadShape]:
        """Enumerate contained solids."""
        from OCP.TopAbs import TopAbs_SOLID
        from OCP.TopExp import TopExp_Explorer

        cast_solid = _topods_cast("Solid")
        out: list[CadShape] = []
        exp = TopExp_Explorer(self.wrapped, TopAbs_SOLID)
        while exp.More():
            out.append(CadShape(cast_solid(exp.Current())))
            exp.Next()
        return out

    # CadQuery compatibility alias — some legacy callers use the camel-cased form.
    def Solids(self) -> list[CadShape]:  # noqa: N802
        """CadQuery-compatible alias for :meth:`solids`."""
        return self.solids()

    def Vertices(self) -> list[tuple[float, float, float]]:  # noqa: N802
        """Enumerate the vertex coordinates of the shape.

        CadQuery-compatible: callers use this to check that a generated mesh
        has any vertices at all (``assert np.any(shape.Vertices())``).
        """
        from OCP.BRep import BRep_Tool
        from OCP.TopAbs import TopAbs_VERTEX
        from OCP.TopExp import TopExp_Explorer

        cast_vertex = _topods_cast("Vertex")
        pnt = getattr(BRep_Tool, "Pnt_s", None) or BRep_Tool.Pnt
        out: list[tuple[float, float, float]] = []
        exp = TopExp_Explorer(self.wrapped, TopAbs_VERTEX)
        while exp.More():
            v = cast_vertex(exp.Current())
            p = pnt(v)
            out.append((float(p.X()), float(p.Y()), float(p.Z())))
            exp.Next()
        return out

    def Faces(self) -> list[CadShape]:  # noqa: N802
        """Enumerate the faces of the shape (CadQuery compatibility)."""
        from OCP.TopAbs import TopAbs_FACE
        from OCP.TopExp import TopExp_Explorer

        cast_face = _topods_cast("Face")
        out: list[CadShape] = []
        exp = TopExp_Explorer(self.wrapped, TopAbs_FACE)
        while exp.More():
            out.append(CadShape(cast_face(exp.Current())))
            exp.Next()
        return out

    def Closed(self) -> bool:  # noqa: N802
        """Whether the shape is topologically closed (CadQuery compatibility).

        Solids are always closed; for shells/compounds we read OCCT's
        per-shape ``Closed`` flag (set by ``BRep_Builder::IsClosed`` when the
        shell was built from a watertight set of faces).
        """
        from OCP.TopAbs import TopAbs_SOLID

        if self.wrapped.ShapeType() == TopAbs_SOLID:
            return True
        return bool(self.wrapped.Closed())

    def Volume(self) -> float:  # noqa: N802
        """Return the (unsigned) volume of the shape.

        OCCT's ``BRepGProp::VolumeProperties`` returns a *signed* volume that
        depends on face orientation; mesh-built shells from
        :func:`mesh_to_shell_brep` can carry inverted orientation and yield a
        negative value.  We return ``abs(...)`` to match CadQuery's behaviour
        and the natural expectation that volumes are non-negative.

        If a mesh-derived volume was stashed on ``_mesh_volume`` AND the OCCT
        solid is not valid (BRepCheck_Analyzer flags self-intersection /
        non-manifold edges, common on raw marching-cubes input), we trust the
        mesh volume — the OCCT surface integral on an invalid topology is
        meaningless.
        """
        from OCP.BRepCheck import BRepCheck_Analyzer
        from OCP.BRepGProp import BRepGProp
        from OCP.GProp import GProp_GProps

        if (
            self._mesh_volume is not None
            and not BRepCheck_Analyzer(
                self.wrapped,
            ).IsValid()
        ):
            return float(abs(self._mesh_volume))

        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(self.wrapped, props)
        return float(abs(props.Mass()))

    def Center(self) -> _Centre:  # noqa: N802
        """Return the volumetric center of mass.

        The result is a :class:`_Centre` — exposes ``.x``, ``.y``, ``.z``,
        ``.toTuple()`` (CadQuery compatibility), and unpacks like a tuple.
        """
        from OCP.BRepGProp import BRepGProp
        from OCP.GProp import GProp_GProps

        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(self.wrapped, props)
        c = props.CentreOfMass()
        return _Centre(float(c.X()), float(c.Y()), float(c.Z()))

    def BoundingBox(self) -> _BBox:  # noqa: N802
        """Return the axis-aligned bounding box.

        The result exposes CadQuery-compatible ``xmin``/``xmax``/… attributes
        (see :class:`_BBox`).
        """
        from OCP.Bnd import Bnd_Box
        from OCP.BRepBndLib import BRepBndLib

        box = Bnd_Box()
        # AddOptimal uses exact geometric bounds (not cached triangulation),
        # matching CadQuery's BoundingBox() behaviour.
        BRepBndLib.AddOptimal_s(self.wrapped, box, True, True)
        xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
        return _BBox(xmin, ymin, zmin, xmax, ymax, zmax)

    # -- exports -----------------------------------------------------------

    def export_stl(
        self,
        path: str | Path,
        linear_deflection: float = 0.01,
        angular_deflection: float = 0.5,
        *,
        ascii_mode: bool = False,
    ) -> None:
        """Export to STL.  Mesh is regenerated at the given deflection."""
        from OCP.BRepMesh import BRepMesh_IncrementalMesh
        from OCP.StlAPI import StlAPI_Writer

        BRepMesh_IncrementalMesh(
            self.wrapped,
            float(linear_deflection),
            False,
            float(angular_deflection),
            True,
        )
        writer = StlAPI_Writer()
        writer.ASCIIMode = bool(ascii_mode)
        if not writer.Write(self.wrapped, str(path)):
            err_msg = f"STL write failed for {path!r}"
            raise RuntimeError(err_msg)

    def export_step(self, path: str | Path) -> None:
        """Export to STEP (AP214)."""
        from OCP.IFSelect import IFSelect_RetDone
        from OCP.STEPControl import (
            STEPControl_AsIs,
            STEPControl_Writer,
        )

        writer = STEPControl_Writer()
        status = writer.Transfer(self.wrapped, STEPControl_AsIs)
        if status != IFSelect_RetDone:
            err_msg = f"STEP transfer failed with status {status!r}"
            raise RuntimeError(err_msg)
        status = writer.Write(str(path))
        if status != IFSelect_RetDone:
            err_msg = f"STEP write failed with status {status!r}"
            raise RuntimeError(err_msg)

    def export_brep(self, path: str | Path) -> None:
        """Export to OCCT native BREP."""
        from OCP.BRepTools import BRepTools

        ok = BRepTools.Write_s(self.wrapped, str(path))
        if not ok:
            err_msg = f"BREP write failed for {path!r}"
            raise RuntimeError(err_msg)


def import_step(path: str | Path) -> CadShape:
    """Import a STEP file and return the resulting :class:`CadShape`.

    Multi-root STEP files are merged into a single ``TopoDS_Compound``.
    """
    require_cad()
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Reader

    reader = STEPControl_Reader()
    status = reader.ReadFile(str(path))
    if status != IFSelect_RetDone:
        err_msg = f"STEP read failed for {path!r} with status {status!r}"
        raise RuntimeError(err_msg)
    reader.TransferRoots()
    return CadShape(reader.OneShape())


# ---------------------------------------------------------------------------
# Mesh → Shell (SOTA path: one TopoDS_Face with attached Poly_Triangulation)
# ---------------------------------------------------------------------------


class ShellCreationError(RuntimeError):
    """Raised when a mesh cannot be converted into an OCCT shell."""


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
    from OCP.BRep import BRep_Builder
    from OCP.gp import gp_Pnt
    from OCP.Poly import Poly_Triangle, Poly_Triangulation
    from OCP.TopoDS import TopoDS_Face, TopoDS_Shell

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


def _triangle_components(
    triangles: npt.NDArray[np.int64],
) -> npt.NDArray[np.int64]:
    """Label connected components of a triangle mesh by edge adjacency."""
    from collections import defaultdict

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
    from collections import defaultdict

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
    from OCP.BRep import BRep_Builder
    from OCP.BRepBuilderAPI import (
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeWire,
    )
    from OCP.gp import gp_Dir, gp_Pln, gp_Pnt
    from OCP.TopoDS import TopoDS_Shell

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
    from OCP.BRep import BRep_Builder
    from OCP.BRepBuilderAPI import (
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeWire,
    )
    from OCP.gp import gp_Pnt
    from OCP.TopoDS import TopoDS_Shell

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
    from OCP.BRepBuilderAPI import (
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeWire,
        BRepBuilderAPI_Sewing,
    )
    from OCP.gp import gp_Pnt

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


# ---------------------------------------------------------------------------
# Primitive builders
# ---------------------------------------------------------------------------


def make_box(dim: Sequence[float], center: Sequence[float]) -> CadShape:
    """Axis-aligned box of size ``dim`` centered at ``center``."""
    require_cad()
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    dx, dy, dz = (float(d) for d in dim)
    cx, cy, cz = (float(c) for c in center)
    corner = gp_Pnt(cx - dx / 2.0, cy - dy / 2.0, cz - dz / 2.0)
    return CadShape(BRepPrimAPI_MakeBox(corner, dx, dy, dz).Shape())


def make_sphere(radius: float, center: Sequence[float]) -> CadShape:
    """Sphere of given radius at ``center``."""
    require_cad()
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere
    from OCP.gp import gp_Pnt

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
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

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
    from OCP.BRepPrimAPI import (
        BRepPrimAPI_MakeCylinder,
        BRepPrimAPI_MakeSphere,
    )
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    cx, cy, cz = (float(c) for c in center)
    h = float(height)
    r = float(radius)
    base_axis = gp_Ax2(gp_Pnt(cx - h / 2.0, cy, cz), gp_Dir(1.0, 0.0, 0.0))
    cyl = BRepPrimAPI_MakeCylinder(base_axis, r, h).Shape()
    left = BRepPrimAPI_MakeSphere(gp_Pnt(cx - h / 2.0, cy, cz), r).Shape()
    right = BRepPrimAPI_MakeSphere(gp_Pnt(cx + h / 2.0, cy, cz), r).Shape()
    fused = CadShape(cyl).fuse(CadShape(left)).fuse(CadShape(right))
    return fused


def make_ellipsoid(radii: Sequence[float], center: Sequence[float]) -> CadShape:
    """Ellipsoid of the given axis-aligned radii at ``center``.

    Built as a unit sphere transformed by a non-uniform scaling matrix.
    """
    require_cad()
    from OCP.BRepBuilderAPI import BRepBuilderAPI_GTransform
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere
    from OCP.gp import gp_GTrsf, gp_Mat, gp_Pnt, gp_XYZ

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
    from OCP.BRepBuilderAPI import (
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeSolid,
        BRepBuilderAPI_MakeWire,
        BRepBuilderAPI_Sewing,
    )
    from OCP.gp import gp_Pnt
    from OCP.ShapeFix import ShapeFix_Solid
    from OCP.TopAbs import TopAbs_SHELL
    from OCP.TopExp import TopExp_Explorer

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
    from OCP.BRepBuilderAPI import (
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeWire,
    )
    from OCP.BRepPrimAPI import BRepPrimAPI_MakePrism
    from OCP.gp import gp_Pnt, gp_Vec

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


# ---------------------------------------------------------------------------
# Compound assembly (for lattice)
# ---------------------------------------------------------------------------


def make_compound(shapes: Iterable[CadShape]) -> CadShape:
    """Assemble shapes into a single OCCT ``TopoDS_Compound``."""
    require_cad()
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound

    builder = BRep_Builder()
    compound = TopoDS_Compound()
    builder.MakeCompound(compound)
    for s in shapes:
        builder.Add(compound, s.wrapped)
    return CadShape(compound)


def make_compound_from_solids(solids: Iterable[Any]) -> CadShape:
    """Assemble raw OCCT ``TopoDS_Shape`` solids (not ``CadShape``) into a compound."""
    require_cad()
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound

    builder = BRep_Builder()
    compound = TopoDS_Compound()
    builder.MakeCompound(compound)
    for s in solids:
        shape = s.wrapped if hasattr(s, "wrapped") else s
        builder.Add(compound, shape)
    return CadShape(compound)


# ---------------------------------------------------------------------------
# Topology exploration and splitting
# ---------------------------------------------------------------------------


def enumerate_solids(shape: CadShape) -> list[Any]:
    """Return the list of ``TopoDS_Solid`` inside a shape (empty if none)."""
    require_cad()
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    cast_solid = _topods_cast("Solid")
    out: list[Any] = []
    exp = TopExp_Explorer(shape.wrapped, TopAbs_SOLID)
    while exp.More():
        out.append(cast_solid(exp.Current()))
        exp.Next()
    return out


def split_shape(
    shape: CadShape,
    tool: CadShape | Iterable[CadShape],
    *,
    fuzzy_value: float = 1e-4,
) -> CadShape:
    """Split *shape* by *tool* using OCCT's ``BRepAlgoAPI_Splitter``.

    The result is a :class:`CadShape` wrapping a ``TopoDS_Compound`` that
    contains the sub-shapes produced by the split.  Use
    :func:`enumerate_solids` to iterate over the resulting solids.

    :param tool: a single :class:`CadShape` or an iterable of them.  Pass
        multiple tools when each tool is a separate ``TopoDS_Shell`` and
        you want to keep them topologically distinct — the OCCT splitter
        treats each list entry as one cutting tool.  This matters because
        ``BRepAlgoAPI_Fuse`` on two shells *decomposes* them into a
        compound of per-face shells, which then defeats the splitter.
    :param fuzzy_value: tolerance forwarded to ``SetFuzzyValue`` so OCCT
        recognises near-coincident geometry as touching.  Necessary when
        a tool is a tessellated shell whose boundary edges lie on
        ``shape``'s planar faces only up to a few microns of drift.
    """
    require_cad()
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Splitter
    from OCP.TopTools import TopTools_ListOfShape

    args = TopTools_ListOfShape()
    args.Append(shape.wrapped)
    tools = TopTools_ListOfShape()
    if isinstance(tool, CadShape):
        tools.Append(tool.wrapped)
    else:
        for t in tool:
            tools.Append(t.wrapped)

    splitter = BRepAlgoAPI_Splitter()
    splitter.SetArguments(args)
    splitter.SetTools(tools)
    if fuzzy_value > 0:
        splitter.SetFuzzyValue(float(fuzzy_value))
    splitter.Build()
    return CadShape(splitter.Shape())


def make_plane_face(
    base_pnt: Sequence[float],
    direction: Sequence[float],
    half_size: float = 1.0e6,
) -> CadShape:
    """Build a large planar face used as a cutting tool.

    :param base_pnt: a point on the plane
    :param direction: plane normal
    :param half_size: half-edge of the square face (default large so the plane
        reaches far outside any realistic shape)
    """
    require_cad()
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    from OCP.gp import gp_Ax3, gp_Dir, gp_Pln, gp_Pnt

    pnt = gp_Pnt(float(base_pnt[0]), float(base_pnt[1]), float(base_pnt[2]))
    nrm = gp_Dir(float(direction[0]), float(direction[1]), float(direction[2]))
    plane = gp_Pln(gp_Ax3(pnt, nrm))
    face = BRepBuilderAPI_MakeFace(
        plane,
        -float(half_size),
        float(half_size),
        -float(half_size),
        float(half_size),
    ).Face()
    return CadShape(face)


def transform_geometry(shape: CadShape, matrix: npt.NDArray[np.float64]) -> CadShape:
    """Apply a 3x4 affine matrix (linear + translation).

    Wraps OCCT ``BRepBuilderAPI_GTransform``.

    :param matrix: ``(3, 4)`` array; rows are ``[a b c tx; d e f ty; g h i tz]``.
    """
    require_cad()
    from OCP.BRepBuilderAPI import BRepBuilderAPI_GTransform
    from OCP.gp import gp_GTrsf, gp_Mat, gp_XYZ

    m = np.asarray(matrix, dtype=np.float64)
    gtrsf = gp_GTrsf()
    gtrsf.SetVectorialPart(
        gp_Mat(
            float(m[0, 0]),
            float(m[0, 1]),
            float(m[0, 2]),
            float(m[1, 0]),
            float(m[1, 1]),
            float(m[1, 2]),
            float(m[2, 0]),
            float(m[2, 1]),
            float(m[2, 2]),
        ),
    )
    gtrsf.SetTranslationPart(gp_XYZ(float(m[0, 3]), float(m[1, 3]), float(m[2, 3])))
    return CadShape(BRepBuilderAPI_GTransform(shape.wrapped, gtrsf, True).Shape())


def translate_solid(solid: Any, offset: Sequence[float]) -> Any:
    """Translate a raw OCCT solid/shape by ``offset`` and return the same type."""
    require_cad()
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.gp import gp_Trsf, gp_Vec

    shape = solid.wrapped if hasattr(solid, "wrapped") else solid
    trsf = gp_Trsf()
    trsf.SetTranslation(
        gp_Vec(float(offset[0]), float(offset[1]), float(offset[2])),
    )
    return BRepBuilderAPI_Transform(shape, trsf, True).Shape()


def solid_center(shape: Any) -> tuple[float, float, float]:
    """Return the volumetric center of mass of a raw OCCT shape/solid."""
    require_cad()
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    s = shape.wrapped if hasattr(shape, "wrapped") else shape
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(s, props)
    com = props.CentreOfMass()
    return (float(com.X()), float(com.Y()), float(com.Z()))


def select_solids_on_side(
    shape: CadShape,
    base_pnt: Sequence[float],
    side_direction: Sequence[float],
) -> list[Any]:
    """Enumerate a compound's solids, keep those whose centroid is on the
    positive side of the plane through ``base_pnt`` normal to ``side_direction``.

    Matches CadQuery's ``.solids(">X")`` / ``.solids("<X")`` semantics for the
    split-by-plane case.  Reverse the sign of ``side_direction`` to get the
    opposite side.
    """
    bx, by, bz = (float(c) for c in base_pnt)
    dx, dy, dz = (float(c) for c in side_direction)
    selected: list[Any] = []
    for solid in enumerate_solids(shape):
        cx, cy, cz = solid_center(solid)
        if (cx - bx) * dx + (cy - by) * dy + (cz - bz) * dz > 0:
            selected.append(solid)
    return selected


def intersect_solids_with_box(solids: Iterable[Any], box: CadShape) -> CadShape:
    """Intersect each solid with *box* and fuse the results into a ``CadShape``.

    Returns a :class:`CadShape` wrapping a compound (possibly empty).
    """
    require_cad()
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common

    parts: list[Any] = []
    for solid in solids:
        s = solid.wrapped if hasattr(solid, "wrapped") else solid
        common = BRepAlgoAPI_Common(s, box.wrapped).Shape()
        # Keep if it contains at least one solid.
        if enumerate_solids(CadShape(common)):
            parts.append(common)
    return make_compound_from_solids(parts)
