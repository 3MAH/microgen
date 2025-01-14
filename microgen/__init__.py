"""Microstructure generation package."""

import importlib.metadata

from .box_mesh import BoxMesh
from .external import Mmg, Neper, parseNeper
from .mesh import is_periodic, mesh, mesh_periodic, meshPeriodic
from .operations import (
    cut_phase_by_shape_list,
    cut_phases,
    cut_phases_by_shape,
    cut_shapes,
    cutPhaseByShapeList,
    cutPhases,
    cutPhasesByShape,
    cutShapes,
    fuse_shapes,
    fuseShapes,
    raster_phase,
    rasterPhase,
    repeat_polydata,
    repeat_shape,
    repeatPolyData,
    repeatShape,
    rescale,
    rotate_euler,
    rotate_pv_euler,
    rotateEuler,
    rotatePvEuler,
)
from .periodic import periodic, periodic_split_and_translate
from .phase import Phase
from .report import Report
from .rve import Rve
from .shape import (
    Box,
    Capsule,
    Cylinder,
    CylindricalTpms,
    Ellipsoid,
    ExtrudedPolygon,
    Infill,
    NormedDistance,
    Polyhedron,
    Sphere,
    SphericalTpms,
    Tpms,
    new_geometry,
    newGeometry,
    surface_functions,
)
from .single_mesh import SingleMesh, check_if_only_linear_tetrahedral

__version__ = importlib.metadata.version(__package__ or __name__)

__all__ = [
    "Box",
    "BoxMesh",
    "Capsule",
    "Cylinder",
    "Ellipsoid",
    "ExtrudedPolygon",
    "Infill",
    "NormedDistance",
    "Mmg",
    "Neper",
    "Phase",
    "Polyhedron",
    "Rve",
    "Sphere",
    "Tpms",
    "CylindricalTpms",
    "SphericalTpms",
    "check_if_only_linear_tetrahedral",
    "cut_phase_by_shape_list",
    "cutPhaseByShapeList",
    "cut_phases",
    "cutPhases",
    "cut_phases_by_shape",
    "cutPhasesByShape",
    "cut_shapes",
    "cutShapes",
    "fuse_shapes",
    "fuseShapes",
    "is_periodic",
    "mesh",
    "mesh_periodic",
    "meshPeriodic",
    "new_geometry",
    "newGeometry",
    "parseNeper",
    "periodic",
    "periodic_split_and_translate",
    "raster_phase",
    "rasterPhase",
    "repeat_polydata",
    "repeatPolyData",
    "repeat_shape",
    "repeatShape",
    "rescale",
    "rotate_euler",
    "rotateEuler",
    "rotate_pv_euler",
    "rotatePvEuler",
    "Report",
    "SingleMesh",
    "surface_functions",
]
