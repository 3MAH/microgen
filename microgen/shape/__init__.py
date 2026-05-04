"""Shape.

========================================
Shape (:mod:`microgen.shape`)
========================================

.. jupyter-execute::
   :hide-code:

   import pyvista
   pyvista.set_jupyter_backend('static')
   pyvista.global_theme.background = 'white'
   pyvista.global_theme.window_size = [600, 400]
   pyvista.global_theme.axes.show = False
   pyvista.global_theme.smooth_shading = True
   pyvista.global_theme.split_sharp_edges = True

"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from . import implicit_ops, surface_functions
from .box import Box
from .capsule import Capsule
from .cylinder import Cylinder
from .ellipsoid import Ellipsoid
from .extruded_polygon import ExtrudedPolygon
from .implicit_ops import batch_smooth_union, from_field
from .polyhedron import Polyhedron
from .shape import Shape
from .sphere import Sphere
from .strut_lattice import (
    AbstractLattice,
    BodyCenteredCubic,
    Cubic,
    Cuboctahedron,
    CustomLattice,
    Diamond,
    FaceCenteredCubic,
    Octahedron,
    OctetTruss,
    RhombicCuboctahedron,
    RhombicDodecahedron,
    TruncatedCube,
    TruncatedCuboctahedron,
    TruncatedOctahedron,
)
from .spinodoid import Spinodoid
from .tpms import CylindricalTpms, Infill, SphericalTpms, Tpms
from .tpms_grading import NormedDistance

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    Vector3DType = tuple[float, float, float] | Sequence[float]

    TpmsPartType = Literal["sheet", "lower skeletal", "upper skeletal", "surface"]

    KwargsGenerateType = int | TpmsPartType

    GeometryParameterType = (
        float
        | Sequence[float]
        | tuple[float, float, float]
        | Sequence[tuple[float, float]]
        | Callable
    )


_SHAPE_FACTORIES: dict[str, Callable[..., Shape]] = {
    "box": lambda center, orientation, p: Box(
        center=center, orientation=orientation, dim=p["dim"],
    ),
    "cylinder": lambda center, orientation, p: Cylinder(
        center=center, orientation=orientation,
        height=p["height"], radius=p["radius"],
    ),
    "extrudedpolygon": lambda center, orientation, p: ExtrudedPolygon(
        center=center, orientation=orientation,
        list_corners=p["list_corners"], height=p["height"],
    ),
    "capsule": lambda center, orientation, p: Capsule(
        center=center, orientation=orientation,
        height=p["height"], radius=p["radius"],
    ),
    "sphere": lambda center, _orientation, p: Sphere(
        center=center, radius=p["radius"],
    ),
    "ellipsoid": lambda center, orientation, p: Ellipsoid(
        center=center, orientation=orientation, radii=p["radii"],
    ),
    "tpms": lambda center, orientation, p: Tpms(
        center=center, orientation=orientation,
        surface_function=p["surface_function"],
        offset=p["offset"], cell_size=p["cell_size"],
        repeat_cell=p["repeat_cell"], resolution=p["resolution"],
    ),
    "polyhedron": lambda _center, _orientation, p: Polyhedron(dic=p["dic"]),
}


def new_geometry(
    shape: str,
    param_geom: dict[str, GeometryParameterType],
    center: tuple[float, float, float] = (0, 0, 0),
    orientation: tuple[float, float, float] = (0, 0, 0),
) -> Shape:
    """Create a new basic geometry with given shape and geometrical parameters.

    :param shape: name of the geometry
    :param param_geom: dictionary with required geometrical parameters
    :param center: center
    :param orientation: orientation

    :return geometry: Shape
    """
    factory = _SHAPE_FACTORIES.get(shape.lower())
    if factory is None:
        raise ShapeError(shape)
    return factory(center, orientation, param_geom)


class ShapeError(Exception):
    """Exception raised for errors in the shape module.

    :param message: explanation of the error
    """

    def __init__(self: ShapeError, shape: str) -> None:
        """Initialize the exception."""
        message = f"{shape} name not implemented"
        super().__init__(message)


__all__ = [
    "AbstractLattice",
    "BodyCenteredCubic",
    "Box",
    "Capsule",
    "Cubic",
    "Cuboctahedron",
    "CustomLattice",
    "Cylinder",
    "CylindricalTpms",
    "Diamond",
    "Ellipsoid",
    "ExtrudedPolygon",
    "FaceCenteredCubic",
    "Infill",
    "NormedDistance",
    "Octahedron",
    "OctetTruss",
    "Polyhedron",
    "RhombicCuboctahedron",
    "RhombicDodecahedron",
    "Shape",
    "Sphere",
    "SphericalTpms",
    "Spinodoid",
    "Tpms",
    "TruncatedCube",
    "TruncatedCuboctahedron",
    "TruncatedOctahedron",
    "batch_smooth_union",
    "from_field",
    "implicit_ops",
    "new_geometry",
    "surface_functions",
]
