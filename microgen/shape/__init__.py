"""
========================================
Shape (:mod:`microgen.shape`)
========================================
"""

from .box import Box
from .capsule import Capsule
from .cylinder import Cylinder
from .ellipsoid import Ellipsoid
from .extrudedPolygon import ExtrudedPolygon
from .polyhedron import Polyhedron
from .sphere import Sphere
from .tpms import Tpms

import cadquery as cq


def newGeometry(
    shape: str,
    center: tuple[float, float, float],
    orientation: tuple[float, float, float],
    param_geom
) -> cq.Shape:
    '''
    Creates a new basic geometry with given shape and geometrical parameters

    :param shape: name of the geometry
    :param center: center
    :param orientation: orientation
    :param param_geom: dictionnary with required geometrical parameters

    :return geometry: cq.Shape
    '''
    if shape.lower() == "box":
        geometry = Box(
            center=center,
            orientation=orientation,
            dim_x=param_geom["dim_x"],
            dim_y=param_geom["dim_y"],
            dim_z=param_geom["dim_z"],
        )
    elif shape.lower() == "cylinder":
        geometry = Cylinder(
            center=center,
            orientation=orientation,
            height=param_geom["height"],
            radius=param_geom["radius"],
        )
    elif shape.lower() == "extrudedpolygon":
        geometry = ExtrudedPolygon(
            center=center,
            orientation=orientation,
            listCorners=param_geom["listCorners"],
            height=param_geom["height"],
        )
    elif shape.lower() == "capsule":
        geometry = Capsule(
            center=center,
            orientation=orientation,
            height=param_geom["height"],
            radius=param_geom["radius"],
        )
    elif shape.lower() == "sphere":
        geometry = Sphere(
            center=center, radius=param_geom["radius"]
        )
    elif shape.lower() == "ellipsoid":
        geometry = Ellipsoid(
            center=center,
            orientation=orientation,
            a_x=param_geom["a_x"],
            a_y=param_geom["a_y"],
            a_z=param_geom["a_z"],
        )
    elif shape.lower() == "tpms":
        geometry = Tpms(
            center=center,
            orientation=orientation,
            surface_function=param_geom["surface_function"],
            type_part=param_geom["type_part"],
            thickness=param_geom["thickness"],
        )
    elif shape.lower() == "polyhedron":
        geometry = Polyhedron(dic=param_geom["dic"])
    else:
        print(shape + ' name not recognised')

    return geometry
