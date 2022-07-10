"""
========================================
Shape (:mod:`microgen.shape`)
========================================

.. jupyter-execute::
   :hide-code:

   import pyvista
   pyvista.set_jupyter_backend('pythreejs')
   pyvista.global_theme.background = 'white'
   pyvista.global_theme.window_size = [690, 400]
   pyvista.global_theme.antialiasing = True
   pyvista.global_theme.axes.show = False

"""

from typing import Any

from .basicGeometry import BasicGeometry
from .box import Box
from .capsule import Capsule
from .cylinder import Cylinder
from .ellipsoid import Ellipsoid
from .extrudedPolygon import ExtrudedPolygon
from .polyhedron import Polyhedron
from .sphere import Sphere
from .tpms import Tpms


def newGeometry(
    shape: str,
    param_geom: dict[str, Any],
    center: tuple[float, float, float] = (0, 0, 0),
    orientation: tuple[float, float, float] = (0, 0, 0),
) -> BasicGeometry:
    """
    Creates a new basic geometry with given shape and geometrical parameters

    :param shape: name of the geometry
    :param param_geom: dictionnary with required geometrical parameters
    :param center: center
    :param orientation: orientation

    :return geometry: BasicGeometry
    """
    if shape.lower() == "box":
        return Box(
            center=center,
            orientation=orientation,
            dim_x=param_geom["dim_x"],
            dim_y=param_geom["dim_y"],
            dim_z=param_geom["dim_z"],
        )
    elif shape.lower() == "cylinder":
        return Cylinder(
            center=center,
            orientation=orientation,
            height=param_geom["height"],
            radius=param_geom["radius"],
        )
    elif shape.lower() == "extrudedpolygon":
        return ExtrudedPolygon(
            center=center,
            orientation=orientation,
            listCorners=param_geom["listCorners"],
            height=param_geom["height"],
        )
    elif shape.lower() == "capsule":
        return Capsule(
            center=center,
            orientation=orientation,
            height=param_geom["height"],
            radius=param_geom["radius"],
        )
    elif shape.lower() == "sphere":
        return Sphere(center=center, radius=param_geom["radius"])
    elif shape.lower() == "ellipsoid":
        return Ellipsoid(
            center=center,
            orientation=orientation,
            a_x=param_geom["a_x"],
            a_y=param_geom["a_y"],
            a_z=param_geom["a_z"],
        )
    elif shape.lower() == "tpms":
        return Tpms(
            center=center,
            orientation=orientation,
            surface_function=param_geom["surface_function"],
            type_part=param_geom["type_part"],
            thickness=param_geom["thickness"],
            cell_size=param_geom["cell_size"],
            repeat_cell=param_geom["repeat_cell"],
        )
    elif shape.lower() == "polyhedron":
        return Polyhedron(dic=param_geom["dic"])
    else:
        raise ValueError(shape + " name not recognised")
