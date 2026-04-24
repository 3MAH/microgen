"""CAD backend — direct OCCT (via OCP) replacement for CadQuery.

=========================================================
CAD backend (:mod:`microgen.cad`)
=========================================================

All CadQuery calls in microgen have been replaced by direct OCP
(``cadquery-ocp``) calls housed in this module.  OCP is an *optional*
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

from typing import TYPE_CHECKING, Iterable, Sequence

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from pathlib import Path

    from OCP.TopoDS import TopoDS_Compound, TopoDS_Shape


_INSTALL_HINT = (
    "microgen's CAD backend requires cadquery-ocp. "
    "Install it with:  pip install 'microgen[cad]'  "
    "or:  pip install cadquery-ocp"
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

    __slots__ = ("xmin", "ymin", "zmin", "xmax", "ymax", "zmax")

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
        import OCP  # noqa: F401, PLC0415
    except ImportError as err:
        raise ImportError(_INSTALL_HINT) from err


# ---------------------------------------------------------------------------
# CadShape wrapper
# ---------------------------------------------------------------------------


class CadShape:
    """Thin wrapper around an OCCT ``TopoDS_Shape``.

    Preserves the ``.wrapped`` attribute name used by CadQuery so downstream
    OCP calls (``BRepAlgoAPI_Fuse(a.wrapped, b.wrapped)``) keep working.
    """

    __slots__ = ("wrapped",)

    def __init__(self, shape: TopoDS_Shape) -> None:
        """Wrap an OCCT ``TopoDS_Shape``."""
        self.wrapped = shape

    # -- transforms --------------------------------------------------------

    def translate(self, offset: Sequence[float]) -> CadShape:
        """Return a translated copy."""
        from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform  # noqa: PLC0415
        from OCP.gp import gp_Trsf, gp_Vec  # noqa: PLC0415

        trsf = gp_Trsf()
        trsf.SetTranslation(gp_Vec(float(offset[0]), float(offset[1]), float(offset[2])))
        transformed = BRepBuilderAPI_Transform(self.wrapped, trsf, True).Shape()
        return CadShape(transformed)

    def rotate(
        self,
        center: Sequence[float],
        axis: Sequence[float],
        angle_degrees: float,
    ) -> CadShape:
        """Return a rotated copy (angle in degrees, axis is a unit vector)."""
        from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform  # noqa: PLC0415
        from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf  # noqa: PLC0415

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
        from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy  # noqa: PLC0415

        return CadShape(BRepBuilderAPI_Copy(self.wrapped).Shape())

    # -- boolean ops -------------------------------------------------------

    def fuse(self, other: CadShape) -> CadShape:
        """Boolean fusion: ``self ∪ other``."""
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse  # noqa: PLC0415

        return CadShape(BRepAlgoAPI_Fuse(self.wrapped, other.wrapped).Shape())

    def cut(self, other: CadShape) -> CadShape:
        """Boolean difference: ``self \\ other``."""
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut  # noqa: PLC0415

        return CadShape(BRepAlgoAPI_Cut(self.wrapped, other.wrapped).Shape())

    def intersect(self, other: CadShape) -> CadShape:
        """Boolean intersection: ``self ∩ other``."""
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Common  # noqa: PLC0415

        return CadShape(BRepAlgoAPI_Common(self.wrapped, other.wrapped).Shape())

    # -- topology queries --------------------------------------------------

    def solids(self) -> list[CadShape]:
        """Enumerate contained solids."""
        from OCP.TopAbs import TopAbs_SOLID  # noqa: PLC0415
        from OCP.TopExp import TopExp_Explorer  # noqa: PLC0415
        from OCP.TopoDS import TopoDS  # noqa: PLC0415

        out: list[CadShape] = []
        exp = TopExp_Explorer(self.wrapped, TopAbs_SOLID)
        while exp.More():
            out.append(CadShape(TopoDS.Solid_s(exp.Current())))
            exp.Next()
        return out

    # CadQuery compatibility alias — some legacy callers use the camel-cased form.
    def Solids(self) -> list[CadShape]:  # noqa: N802
        """CadQuery-compatible alias for :meth:`solids`."""
        return self.solids()

    def Volume(self) -> float:  # noqa: N802
        """Return the volume of the shape (uses OCCT ``BRepGProp``)."""
        from OCP.BRepGProp import BRepGProp  # noqa: PLC0415
        from OCP.GProp import GProp_GProps  # noqa: PLC0415

        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(self.wrapped, props)
        return float(props.Mass())

    def Center(self) -> _Centre:  # noqa: N802
        """Return the volumetric center of mass.

        The result is a :class:`_Centre` — exposes ``.x``, ``.y``, ``.z``,
        ``.toTuple()`` (CadQuery compatibility), and unpacks like a tuple.
        """
        from OCP.BRepGProp import BRepGProp  # noqa: PLC0415
        from OCP.GProp import GProp_GProps  # noqa: PLC0415

        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(self.wrapped, props)
        c = props.CentreOfMass()
        return _Centre(float(c.X()), float(c.Y()), float(c.Z()))

    def BoundingBox(self) -> _BBox:  # noqa: N802
        """Return the axis-aligned bounding box.

        The result exposes CadQuery-compatible ``xmin``/``xmax``/… attributes
        (see :class:`_BBox`).
        """
        from OCP.Bnd import Bnd_Box  # noqa: PLC0415
        from OCP.BRepBndLib import BRepBndLib  # noqa: PLC0415

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
        from OCP.BRepMesh import BRepMesh_IncrementalMesh  # noqa: PLC0415
        from OCP.StlAPI import StlAPI_Writer  # noqa: PLC0415

        BRepMesh_IncrementalMesh(
            self.wrapped,
            float(linear_deflection),
            False,
            float(angular_deflection),
            True,
        )
        writer = StlAPI_Writer()
        writer.ASCIIMode = bool(ascii_mode)
        writer.Write(self.wrapped, str(path))

    def export_step(self, path: str | Path) -> None:
        """Export to STEP (AP214)."""
        from OCP.IFSelect import IFSelect_RetDone  # noqa: PLC0415
        from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer  # noqa: PLC0415

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
        from OCP.BRepTools import BRepTools  # noqa: PLC0415

        ok = BRepTools.Write_s(self.wrapped, str(path))
        if not ok:
            err_msg = f"BREP write failed for {path!r}"
            raise RuntimeError(err_msg)


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
        triangulation.SetNode(i + 1, gp_Pnt(float(pts[i, 0]), float(pts[i, 1]), float(pts[i, 2])))
    for i in range(nb_tri):
        a, b, c = int(tris[i, 0]), int(tris[i, 1]), int(tris[i, 2])
        triangulation.SetTriangle(i + 1, Poly_Triangle(a + 1, b + 1, c + 1))

    builder = BRep_Builder()
    face = TopoDS_Face()
    try:
        builder.MakeFace(face, triangulation)
    except Exception as err:  # noqa: BLE001
        err_msg = "OCCT refused the triangulation — check bounds and field."
        raise ShellCreationError(err_msg) from err

    shell = TopoDS_Shell()
    builder.MakeShell(shell)
    builder.Add(shell, face)
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
    except Exception as err:  # noqa: BLE001
        err_msg = (
            "Failed to build the OCCT shell from the mesh; "
            "try to increase the resolution or adjust bounds."
        )
        raise ShellCreationError(err_msg) from err

    return CadShape(shell)


# ---------------------------------------------------------------------------
# Primitive builders
# ---------------------------------------------------------------------------


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
    fused = CadShape(cyl).fuse(CadShape(left)).fuse(CadShape(right))
    return fused


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
            rx, 0.0, 0.0,
            0.0, ry, 0.0,
            0.0, 0.0, rz,
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
    from OCP.TopoDS import TopoDS  # noqa: PLC0415

    cx, cy, cz = (float(c) for c in center)
    points = [
        gp_Pnt(float(v[0]) + cx, float(v[1]) + cy, float(v[2]) + cz)
        for v in vertices
    ]

    sewing = BRepBuilderAPI_Sewing()
    for ixs in faces_ixs:
        wire_builder = BRepBuilderAPI_MakeWire()
        for i1, i2 in zip(ixs, ixs[1:]):
            edge = BRepBuilderAPI_MakeEdge(points[i1], points[i2]).Edge()
            wire_builder.Add(edge)
        face = BRepBuilderAPI_MakeFace(wire_builder.Wire()).Face()
        sewing.Add(face)
    sewing.Perform()
    sewn = sewing.SewedShape()

    # Extract the shell (sewing may return it directly or wrapped in a compound).
    exp = TopExp_Explorer(sewn, TopAbs_SHELL)
    if not exp.More():
        err_msg = "Sewing did not produce a shell — check face connectivity"
        raise ShellCreationError(err_msg)
    shell = TopoDS.Shell_s(exp.Current())

    solid = BRepBuilderAPI_MakeSolid(shell).Solid()
    fixer = ShapeFix_Solid(solid)
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
        edge = BRepBuilderAPI_MakeEdge(pts[i], pts[i + 1]).Edge()
        wire_builder.Add(edge)
    face = BRepBuilderAPI_MakeFace(wire_builder.Wire()).Face()
    extruded = BRepPrimAPI_MakePrism(face, gp_Vec(h, 0.0, 0.0)).Shape()
    return CadShape(extruded)


# ---------------------------------------------------------------------------
# Compound assembly (for lattice)
# ---------------------------------------------------------------------------


def make_compound(shapes: Iterable[CadShape]) -> CadShape:
    """Assemble shapes into a single OCCT ``TopoDS_Compound``."""
    require_cad()
    from OCP.BRep import BRep_Builder  # noqa: PLC0415
    from OCP.TopoDS import TopoDS_Compound  # noqa: PLC0415

    builder = BRep_Builder()
    compound = TopoDS_Compound()
    builder.MakeCompound(compound)
    for s in shapes:
        builder.Add(compound, s.wrapped)
    return CadShape(compound)


def make_compound_from_solids(solids: Iterable[Any]) -> CadShape:
    """Assemble raw OCCT ``TopoDS_Shape`` solids (not ``CadShape``) into a compound."""
    require_cad()
    from OCP.BRep import BRep_Builder  # noqa: PLC0415
    from OCP.TopoDS import TopoDS_Compound  # noqa: PLC0415

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
    from OCP.TopAbs import TopAbs_SOLID  # noqa: PLC0415
    from OCP.TopExp import TopExp_Explorer  # noqa: PLC0415
    from OCP.TopoDS import TopoDS  # noqa: PLC0415

    out: list[Any] = []
    exp = TopExp_Explorer(shape.wrapped, TopAbs_SOLID)
    while exp.More():
        out.append(TopoDS.Solid_s(exp.Current()))
        exp.Next()
    return out


def split_shape(shape: CadShape, tool: CadShape) -> CadShape:
    """Split *shape* by *tool* using OCCT's ``BRepAlgoAPI_Splitter``.

    The result is a :class:`CadShape` wrapping a ``TopoDS_Compound`` that
    contains the sub-shapes produced by the split.  Use
    :func:`enumerate_solids` to iterate over the resulting solids.
    """
    require_cad()
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Splitter  # noqa: PLC0415
    from OCP.TopTools import TopTools_ListOfShape  # noqa: PLC0415

    args = TopTools_ListOfShape()
    args.Append(shape.wrapped)
    tools = TopTools_ListOfShape()
    tools.Append(tool.wrapped)

    splitter = BRepAlgoAPI_Splitter()
    splitter.SetArguments(args)
    splitter.SetTools(tools)
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
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace  # noqa: PLC0415
    from OCP.gp import gp_Ax3, gp_Dir, gp_Pln, gp_Pnt  # noqa: PLC0415

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
    """Apply a 3x4 affine matrix (linear + translation) via ``BRepBuilderAPI_GTransform``.

    :param matrix: ``(3, 4)`` array; rows are ``[a b c tx; d e f ty; g h i tz]``.
    """
    require_cad()
    from OCP.BRepBuilderAPI import BRepBuilderAPI_GTransform  # noqa: PLC0415
    from OCP.gp import gp_GTrsf, gp_Mat, gp_XYZ  # noqa: PLC0415

    m = np.asarray(matrix, dtype=np.float64)
    gtrsf = gp_GTrsf()
    gtrsf.SetVectorialPart(
        gp_Mat(
            float(m[0, 0]), float(m[0, 1]), float(m[0, 2]),
            float(m[1, 0]), float(m[1, 1]), float(m[1, 2]),
            float(m[2, 0]), float(m[2, 1]), float(m[2, 2]),
        ),
    )
    gtrsf.SetTranslationPart(gp_XYZ(float(m[0, 3]), float(m[1, 3]), float(m[2, 3])))
    return CadShape(BRepBuilderAPI_GTransform(shape.wrapped, gtrsf, True).Shape())


def translate_solid(solid: Any, offset: Sequence[float]) -> Any:
    """Translate a raw OCCT solid/shape by ``offset`` and return the same type."""
    require_cad()
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform  # noqa: PLC0415
    from OCP.gp import gp_Trsf, gp_Vec  # noqa: PLC0415

    shape = solid.wrapped if hasattr(solid, "wrapped") else solid
    trsf = gp_Trsf()
    trsf.SetTranslation(
        gp_Vec(float(offset[0]), float(offset[1]), float(offset[2])),
    )
    return BRepBuilderAPI_Transform(shape, trsf, True).Shape()


def solid_center(shape: Any) -> tuple[float, float, float]:
    """Return the volumetric center of mass of a raw OCCT shape/solid."""
    require_cad()
    from OCP.BRepGProp import BRepGProp  # noqa: PLC0415
    from OCP.GProp import GProp_GProps  # noqa: PLC0415

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
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common  # noqa: PLC0415

    parts: list[Any] = []
    for solid in solids:
        s = solid.wrapped if hasattr(solid, "wrapped") else solid
        common = BRepAlgoAPI_Common(s, box.wrapped).Shape()
        # Keep if it contains at least one solid.
        if enumerate_solids(CadShape(common)):
            parts.append(common)
    return make_compound_from_solids(parts)
