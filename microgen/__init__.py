import importlib.metadata

from .box_mesh import BoxMesh
from .external import Mmg, Neper, parseNeper
from .mesh import is_periodic, mesh, meshPeriodic
from .operations import (
    cutPhaseByShapeList,
    cutPhases,
    cutPhasesByShape,
    cutShapes,
    fuseShapes,
    rasterPhase,
    repeatPolyData,
    repeatShape,
    rescale,
    rotateEuler,
    rotatePvEuler,
)
from .periodic import periodic
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
    Polyhedron,
    Sphere,
    SphericalTpms,
    Tpms,
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
    "cutPhaseByShapeList",
    "cutPhases",
    "cutPhasesByShape",
    "cutShapes",
    "fuseShapes",
    "is_periodic",
    "mesh",
    "meshPeriodic",
    "newGeometry",
    "parseNeper",
    "periodic",
    "rasterPhase",
    "repeatPolyData",
    "repeatShape",
    "rescale",
    "rotateEuler",
    "rotatePvEuler",
    "Report",
    "SingleMesh",
    "surface_functions",
]
