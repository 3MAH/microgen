"""Topology assembly, exploration, and splitting helpers.

Free functions for working with OCCT ``TopoDS_Solid`` / ``TopoDS_Compound``
collections: building compounds, enumerating solids, splitting shapes by
tool, plane-face cutting tools, affine transforms, side selection, and
box-clipping a list of solids.

These power the periodic split-and-translate algorithm in
``microgen.periodic`` and the rasterisation pipeline on ``Phase``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
import numpy.typing as npt

from ._install import require_cad
from .shape import CadShape, _topods_cast

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence


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


def enumerate_solids(shape: CadShape) -> list[Any]:
    """Return the list of ``TopoDS_Solid`` inside a shape (empty if none)."""
    require_cad()
    from OCP.TopAbs import TopAbs_SOLID  # noqa: PLC0415
    from OCP.TopExp import TopExp_Explorer  # noqa: PLC0415

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
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Splitter  # noqa: PLC0415
    from OCP.TopTools import TopTools_ListOfShape  # noqa: PLC0415

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
    """Apply a 3x4 affine matrix (linear + translation).

    Wraps OCCT ``BRepBuilderAPI_GTransform``.

    :param matrix: ``(3, 4)`` array; rows are ``[a b c tx; d e f ty; g h i tz]``.
    """
    require_cad()
    from OCP.BRepBuilderAPI import BRepBuilderAPI_GTransform  # noqa: PLC0415
    from OCP.gp import gp_GTrsf, gp_Mat, gp_XYZ  # noqa: PLC0415

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
    """Enumerate a compound's solids; keep those on the positive side of a plane.

    The plane passes through ``base_pnt`` normal to ``side_direction``.
    Matches CadQuery's ``.solids(">X")`` / ``.solids("<X")`` semantics for
    the split-by-plane case.  Reverse the sign of ``side_direction`` to
    get the opposite side.
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
