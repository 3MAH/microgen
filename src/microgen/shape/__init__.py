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

from typing import TYPE_CHECKING, Callable, Literal, Sequence, Tuple

from . import surface_functions
from .box import Box
from .capsule import Capsule
from .cylinder import Cylinder
from .ellipsoid import Ellipsoid
from .extruded_polygon import ExtrudedPolygon
from .polyhedron import Polyhedron
from .shape import Shape
from .sphere import Sphere
from .tpms import CylindricalTpms, Infill, SphericalTpms, Tpms
from .tpms_grading import NormedDistance

if TYPE_CHECKING:
    Vector3DType = Tuple[float, float, float] | Sequence[float]

    TpmsPartType = Literal["sheet", "lower skeletal", "upper skeletal", "surface"]

    KwargsGenerateType = int | TpmsPartType

    GeometryParameterType = (
        float
        | Sequence[float]
        | tuple[float, float, float]
        | Sequence[tuple[float, float]]
        | Callable
    )


def new_geometry(  # noqa: PLR0911
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
    if shape.lower() == "box":
        return Box(
            center=center,
            orientation=orientation,
            dim=param_geom["dim"],
        )
    if shape.lower() == "cylinder":
        return Cylinder(
            center=center,
            orientation=orientation,
            height=param_geom["height"],
            radius=param_geom["radius"],
        )
    if shape.lower() == "extrudedpolygon":
        return ExtrudedPolygon(
            center=center,
            orientation=orientation,
            listCorners=param_geom["listCorners"],
            height=param_geom["height"],
        )
    if shape.lower() == "capsule":
        return Capsule(
            center=center,
            orientation=orientation,
            height=param_geom["height"],
            radius=param_geom["radius"],
        )
    if shape.lower() == "sphere":
        return Sphere(center=center, radius=param_geom["radius"])
    if shape.lower() == "ellipsoid":
        return Ellipsoid(
            center=center,
            orientation=orientation,
            radii=param_geom["radii"],
        )
    if shape.lower() == "tpms":
        return Tpms(
            center=center,
            orientation=orientation,
            surface_function=param_geom["surface_function"],
            offset=param_geom["offset"],
            cell_size=param_geom["cell_size"],
            repeat_cell=param_geom["repeat_cell"],
            resolution=param_geom["resolution"],
        )
    if shape.lower() == "polyhedron":
        return Polyhedron(dic=param_geom["dic"])

    raise ShapeError(shape)


class ShapeError(Exception):
    """Exception raised for errors in the shape module.

    :param message: explanation of the error
    """

    def __init__(self: ShapeError, shape: str) -> None:
        """Initialize the exception."""
        message = f"{shape} name not implemented"
        super().__init__(message)


# Deprecated
newGeometry = new_geometry  # noqa: N816

__all__ = [
    "Box",
    "Capsule",
    "CylindricalTpms",
    "Cylinder",
    "Ellipsoid",
    "ExtrudedPolygon",
    "Infill",
    "NormedDistance",
    "Polyhedron",
    "Shape",
    "SphericalTpms",
    "Sphere",
    "Tpms",
    "new_geometry",
    "newGeometry",
    "surface_functions",
]
