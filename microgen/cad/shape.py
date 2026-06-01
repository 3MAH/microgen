"""``CadShape`` wrapper and OCCT topology helpers.

This module holds the ``CadShape`` thin wrapper around ``TopoDS_Shape`` plus
the low-level OCP utilities used by the rest of the CAD subpackage
(``_run_boolean``, ``_topods_cast``, ``_Centre``, ``_BBox``,
``ShellCreationError``).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from ._install import require_cad  # noqa: F401  (re-exported for back-compat)

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path

    from OCP.TopoDS import TopoDS_Shape


class _Centre(tuple):
    """Tuple-like 3D point exposing ``.x``, ``.y``, ``.z`` and ``.to_tuple()``.

    Returned by :meth:`CadShape.center`.  Mimics just enough of
    ``cadquery.Vector`` to drop-in where existing code uses
    ``shape.center().to_tuple()`` or ``shape.center().x``.
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

    def to_tuple(self) -> tuple[float, float, float]:
        """Return ``(x, y, z)`` as a plain tuple."""
        return (self[0], self[1], self[2])


class _BBox:
    """Axis-aligned bounding box exposing ``xmin`` / ``xmax`` / ….

    Returned by :meth:`CadShape.bounding_box`.  Also indexable as a 6-tuple
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
    def diagonal_length(self) -> float:
        """Length of the box's space diagonal."""
        dx = self.xmax - self.xmin
        dy = self.ymax - self.ymin
        dz = self.zmax - self.zmin
        return float((dx * dx + dy * dy + dz * dz) ** 0.5)


class ShellCreationError(RuntimeError):
    """Raised when a mesh cannot be converted into an OCCT shell."""


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
    from OCP.TopoDS import TopoDS  # noqa: PLC0415

    return getattr(TopoDS, f"{name}_s", None) or getattr(TopoDS, name)


class CadShape:
    """Thin wrapper around an OCCT ``TopoDS_Shape``.

    Preserves the ``.wrapped`` attribute name used by CadQuery so downstream
    OCP calls (``BRepAlgoAPI_Fuse(a.wrapped, b.wrapped)``) keep working.

    ``_mesh_volume`` (optional) is a trusted volume in the source mesh's
    units, set by mesh-derived constructors (e.g. the TPMS periodic shell)
    where OCCT's surface-integral volume is unreliable on invalid topology.
    :meth:`volume` prefers it over the OCCT integral when present.
    """

    __slots__ = ("_mesh_volume", "wrapped")

    def __init__(self, shape: TopoDS_Shape) -> None:
        """Wrap an OCCT ``TopoDS_Shape``."""
        self.wrapped = shape
        self._mesh_volume: float | None = None

    # -- transforms --------------------------------------------------------

    def translate(self, offset: Sequence[float]) -> CadShape:
        """Return a translated copy."""
        from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform  # noqa: PLC0415
        from OCP.gp import gp_Trsf, gp_Vec  # noqa: PLC0415

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

        return CadShape(_run_boolean(BRepAlgoAPI_Fuse, self, other, "Fuse"))

    def cut(self, other: CadShape) -> CadShape:
        """Boolean difference: ``self \\ other``."""
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut  # noqa: PLC0415

        return CadShape(_run_boolean(BRepAlgoAPI_Cut, self, other, "Cut"))

    def intersect(self, other: CadShape) -> CadShape:
        """Boolean intersection: ``self ∩ other``."""
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Common  # noqa: PLC0415

        return CadShape(_run_boolean(BRepAlgoAPI_Common, self, other, "Common"))

    # -- topology queries --------------------------------------------------

    def solids(self) -> list[CadShape]:
        """Enumerate contained solids."""
        from OCP.TopAbs import TopAbs_SOLID  # noqa: PLC0415
        from OCP.TopExp import TopExp_Explorer  # noqa: PLC0415

        cast_solid = _topods_cast("Solid")
        out: list[CadShape] = []
        exp = TopExp_Explorer(self.wrapped, TopAbs_SOLID)
        while exp.More():
            out.append(CadShape(cast_solid(exp.Current())))
            exp.Next()
        return out

    def vertices(self) -> list[tuple[float, float, float]]:
        """Enumerate the vertex coordinates of the shape.

        Callers use this to check that a generated mesh has any vertices at
        all (``assert np.any(shape.vertices())``).
        """
        from OCP.BRep import BRep_Tool  # noqa: PLC0415
        from OCP.TopAbs import TopAbs_VERTEX  # noqa: PLC0415
        from OCP.TopExp import TopExp_Explorer  # noqa: PLC0415

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

    def faces(self) -> list[CadShape]:
        """Enumerate the faces of the shape."""
        from OCP.TopAbs import TopAbs_FACE  # noqa: PLC0415
        from OCP.TopExp import TopExp_Explorer  # noqa: PLC0415

        cast_face = _topods_cast("Face")
        out: list[CadShape] = []
        exp = TopExp_Explorer(self.wrapped, TopAbs_FACE)
        while exp.More():
            out.append(CadShape(cast_face(exp.Current())))
            exp.Next()
        return out

    def is_closed(self) -> bool:
        """Whether the shape is topologically closed.

        Solids are always closed; for shells/compounds we read OCCT's
        per-shape ``Closed`` flag (set by ``BRep_Builder::IsClosed`` when the
        shell was built from a watertight set of faces).
        """
        from OCP.TopAbs import TopAbs_SOLID  # noqa: PLC0415

        if self.wrapped.ShapeType() == TopAbs_SOLID:
            return True
        return bool(self.wrapped.Closed())

    def volume(self) -> float:
        """Return the (unsigned) volume of the shape.

        OCCT's ``BRepGProp::VolumeProperties`` returns a *signed* volume that
        depends on face orientation; mesh-built shells from
        :func:`microgen.cad.mesh_to_shell_brep` can carry inverted orientation
        and yield a negative value.  We return ``abs(...)`` so volumes are
        non-negative.

        If a mesh-derived volume was stashed on ``_mesh_volume`` AND the OCCT
        solid is not valid (BRepCheck_Analyzer flags self-intersection /
        non-manifold edges, common on raw marching-cubes input), we trust the
        mesh volume — the OCCT surface integral on an invalid topology is
        meaningless.
        """
        from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: PLC0415
        from OCP.BRepGProp import BRepGProp  # noqa: PLC0415
        from OCP.GProp import GProp_GProps  # noqa: PLC0415

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

    def center(self) -> _Centre:
        """Return the volumetric center of mass.

        The result is a :class:`_Centre` — exposes ``.x``, ``.y``, ``.z``,
        ``.to_tuple()``, and unpacks like a tuple.
        """
        from OCP.BRepGProp import BRepGProp  # noqa: PLC0415
        from OCP.GProp import GProp_GProps  # noqa: PLC0415

        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(self.wrapped, props)
        c = props.CentreOfMass()
        return _Centre(float(c.X()), float(c.Y()), float(c.Z()))

    def bounding_box(self) -> _BBox:
        """Return the axis-aligned bounding box.

        The result exposes ``xmin`` / ``xmax`` / … attributes (see
        :class:`_BBox`).
        """
        from OCP.Bnd import Bnd_Box  # noqa: PLC0415
        from OCP.BRepBndLib import BRepBndLib  # noqa: PLC0415

        box = Bnd_Box()
        # AddOptimal uses exact geometric bounds (not cached triangulation).
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
        if not writer.Write(self.wrapped, str(path)):
            err_msg = f"STL write failed for {path!r}"
            raise RuntimeError(err_msg)

    def export_step(self, path: str | Path) -> None:
        """Export to STEP (AP214)."""
        from OCP.IFSelect import IFSelect_RetDone  # noqa: PLC0415
        from OCP.STEPControl import (  # noqa: PLC0415
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
        from OCP.BRepTools import BRepTools  # noqa: PLC0415

        ok = BRepTools.Write_s(self.wrapped, str(path))
        if not ok:
            err_msg = f"BREP write failed for {path!r}"
            raise RuntimeError(err_msg)
