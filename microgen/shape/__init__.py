"""
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

from typing import Any, Dict, Tuple

from . import surface_functions
from .basicGeometry import BasicGeometry
from .box import Box
from .capsule import Capsule
from .cylinder import Cylinder
from .ellipsoid import Ellipsoid
from .extrudedPolygon import ExtrudedPolygon
from .polyhedron import Polyhedron
from .sphere import Sphere
from .tpms import CylindricalTpms, SphericalTpms, Tpms


def newGeometry(
    shape: str,
    param_geom: Dict[str, Any],
    center: Tuple[float, float, float] = (0, 0, 0),
    orientation: Tuple[float, float, float] = (0, 0, 0),
) -> BasicGeometry:
    """
    Creates a new basic geometry with given shape and geometrical parameters

    :param shape: name of the geometry
    :param param_geom: dictionary with required geometrical parameters
    :param center: center
    :param orientation: orientation

    :return geometry: BasicGeometry
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

    raise ValueError(f"{shape} name not recognised")


__all__ = [
    "Box",
    "Capsule",
    "CylindricalTpms",
    "Cylinder",
    "Ellipsoid",
    "ExtrudedPolygon",
    "Polyhedron",
    "SphericalTpms",
    "Sphere",
    "Tpms",
    "newGeometry",
    "surface_functions",
]
