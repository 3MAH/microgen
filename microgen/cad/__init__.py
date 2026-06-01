"""CAD backend — direct OCCT (via OCP) replacement for CadQuery.

=========================================================
CAD backend (:mod:`microgen.cad`)
=========================================================

All CadQuery calls in microgen have been replaced by direct OCP
(``cadquery-ocp-novtk``) calls housed in this subpackage.  OCP is an
*optional* dependency — install via ``pip install microgen[cad]``.

The subpackage's submodules don't import OCP at top level, so
``import microgen.cad`` always succeeds.  OCP-dependent functions import
it lazily and raise a helpful ``ImportError`` with install instructions if
OCP is missing.

Return type
-----------

CAD-producing functions return a :class:`CadShape` — a thin wrapper around
an OCCT ``TopoDS_Shape`` exposing ``.wrapped`` for downstream OCP calls, plus
convenience methods (``translate``, ``rotate``, ``fuse``, ``cut``,
``export_stl``, ``export_step``, ``export_brep``).  ``.wrapped`` matches the
attribute name CadQuery's ``Shape`` exposed, so most legacy call sites keep
working unchanged.

Layout
------

- ``cad._install`` — ``require_cad`` / ``_INSTALL_HINT`` gate
- ``cad.shape``    — ``CadShape`` class + ``_Centre`` / ``_BBox`` /
  ``_run_boolean`` / ``_topods_cast`` / ``ShellCreationError``
- ``cad.io``       — ``import_step`` (export methods live on ``CadShape``)
- ``cad.meshbridge`` — ``mesh_to_shape`` / ``shape_to_cad`` /
  ``mesh_to_planar_face`` / ``mesh_to_shell_brep`` / ``mesh_to_sewn_shell``
- ``cad.primitives`` — ``make_box`` / ``make_sphere`` / ``make_cylinder`` /
  ``make_capsule`` / ``make_ellipsoid`` / ``make_polyhedron`` /
  ``make_extruded_polygon``
- ``cad.topo``     — ``make_compound`` / ``make_compound_from_solids`` /
  ``enumerate_solids`` / ``split_shape`` / ``make_plane_face`` /
  ``transform_geometry`` / ``translate_solid`` / ``solid_center`` /
  ``select_solids_on_side`` / ``intersect_solids_with_box``

All public symbols are re-exported here so existing
``from microgen.cad import …`` imports keep working unchanged.
"""

from __future__ import annotations

from ._install import _INSTALL_HINT, require_cad
from .io import import_step
from .meshbridge import (
    _triangle_components,
    _walk_boundary_loops,
    mesh_to_planar_face,
    mesh_to_sewn_shell,
    mesh_to_shape,
    mesh_to_shell_brep,
    shape_to_cad,
)
from .primitives import (
    make_box,
    make_capsule,
    make_cylinder,
    make_ellipsoid,
    make_extruded_polygon,
    make_polyhedron,
    make_sphere,
)
from .shape import (
    CadShape,
    ShellCreationError,
    _BBox,
    _Centre,
    _run_boolean,
    _topods_cast,
)
from .topo import (
    enumerate_solids,
    intersect_solids_with_box,
    make_compound,
    make_compound_from_solids,
    make_plane_face,
    select_solids_on_side,
    solid_center,
    split_shape,
    transform_geometry,
    translate_solid,
)

__all__ = [
    "CadShape",
    "ShellCreationError",
    "_BBox",
    "_Centre",
    "_INSTALL_HINT",
    "_run_boolean",
    "_topods_cast",
    "_triangle_components",
    "_walk_boundary_loops",
    "enumerate_solids",
    "import_step",
    "intersect_solids_with_box",
    "make_box",
    "make_capsule",
    "make_compound",
    "make_compound_from_solids",
    "make_cylinder",
    "make_ellipsoid",
    "make_extruded_polygon",
    "make_plane_face",
    "make_polyhedron",
    "make_sphere",
    "mesh_to_planar_face",
    "mesh_to_sewn_shell",
    "mesh_to_shape",
    "mesh_to_shell_brep",
    "require_cad",
    "select_solids_on_side",
    "shape_to_cad",
    "solid_center",
    "split_shape",
    "transform_geometry",
    "translate_solid",
]
