"""
====================================================
Basic Geometry (:mod:`microgen.shape.basicGeometry`)
====================================================
"""

import cadquery as cq

from typing import Union

# import box
# from .capsule import Capsule
# from .cylinder import Cylinder
# from .ellipsoid import Ellipsoid
# from .extrudedPolygon import ExtrudedPolygon
# from .polyhedron import Polyhedron
# from .sphere import Sphere
# from .tpms import Tpms


class BasicGeometry:
    """
    BasicGeometry class to manage shapes

    :param shape: name of the shape
    :param center: center
    :param orientation: orientation
    """

    numInstances = 0

    def __init__(
        self,
        shape: str,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
    ) -> None:

        self.number = self.numInstances
        self.shape = shape
        self.center = center
        self.orientation = orientation
        self.name = self.shape + "_" + str(self.number)

        self.geometry = None  # type: Union[cq.Shape, None]
        BasicGeometry.numInstances += 1
