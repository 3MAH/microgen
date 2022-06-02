"""
BasicGeometry class to manage shapes
"""
from typing import Any

import cadquery as cq
import numpy as np

from .rve import Rve
from .shape.box import Box
from .shape.capsule import Capsule
from .shape.cylinder import Cylinder
from .shape.ellipsoid import Ellipsoid
from .shape.extrudedPolygon import ExtrudedPolygon
from .shape.polyhedron import Polyhedron
from .shape.sphere import Sphere
from .shape.tpms import Tpms


class BasicGeometry:
    """
        BasicGeometry class to manage shapes

        :param shape: name of the shape
        :type shape: str
        :param param_geom: dictionary containing corresponding shape parameters
        :type param_geom: dict
    """
    def __init__(
        self,
        shape: str,
        param_geom: dict,
        number: int = 0,
        xc: float = 0,
        yc: float = 0,
        zc: float = 0,
        psi: float = 0,
        theta: float = 0,
        phi: float = 0,
        path_data: str = None,
    ) -> None:
        self.number = number
        self.shape = shape
        self.xc = xc
        self.yc = yc
        self.zc = zc
        self.psi = psi
        self.theta = theta
        self.phi = phi
        self.param_geom = param_geom
        self.path_data = path_data

        self.center = np.array([self.xc, self.yc, self.zc])
        self.orientation = np.array([self.psi, self.theta, self.phi])
        self.name = self.shape + str(self.number)

        self.geometry = None  # type: Any

        if self.shape.lower() == "box":
            self.geometry = Box(
                center=self.center,
                orientation=self.orientation,
                dim_x=self.param_geom["dim_x"],
                dim_y=self.param_geom["dim_y"],
                dim_z=self.param_geom["dim_z"],
                number=self.number,
            )
        if self.shape.lower() == "cylinder":
            self.geometry = Cylinder(
                center=self.center,
                orientation=self.orientation,
                height=self.param_geom["height"],
                radius=self.param_geom["radius"],
                number=self.number,
            )
        if self.shape.lower() == "extrudedpolygon":
            self.geometry = ExtrudedPolygon(
                center=self.center,
                orientation=self.orientation,
                listCorners=self.param_geom["listCorners"],
                height=self.param_geom["height"],
                number=self.number,
            )
        if self.shape.lower() == "capsule":
            self.geometry = Capsule(
                center=self.center,
                orientation=self.orientation,
                height=self.param_geom["height"],
                radius=self.param_geom["radius"],
                number=self.number,
            )
        if self.shape.lower() == "sphere":
            self.geometry = Sphere(
                center=self.center, radius=self.param_geom["radius"], number=self.number
            )
        if self.shape.lower() == "ellipsoid":
            self.geometry = Ellipsoid(
                center=self.center,
                orientation=self.orientation,
                a_x=self.param_geom["a_x"],
                a_y=self.param_geom["a_y"],
                a_z=self.param_geom["a_z"],
                number=self.number,
            )
        if self.shape.lower() == "tpms":
            self.geometry = Tpms(
                center=self.center,
                orientation=self.orientation,
                surface_function=self.param_geom["surface_function"],
                type_part=self.param_geom["type_part"],
                thickness=self.param_geom["thickness"],
                number=self.number,
            )
            # if self.param_geom["type_surface"] == "custom":
            #     self.geometry = Tpms(
            #         center=self.center,
            #         orientation=self.orientation,
            #         type_surface=self.param_geom["type_surface"],
            #         type_part=self.param_geom["type_part"],
            #         thickness=self.param_geom["thickness"],
            #         number=self.number,
            #         function=self.param_geom["function"],
            #     )
            # else:
            #     self.geometry = Tpms(
            #         center=self.center,
            #         orientation=self.orientation,
            #         type_surface=self.param_geom["type_surface"],
            #         type_part=self.param_geom["type_part"],
            #         thickness=self.param_geom["thickness"],
            #         number=self.number,
            #     )

        if self.shape.lower() == "polyhedron":
            self.geometry = Polyhedron(dic=self.param_geom["dic"], number=self.number)

    def __cmp__(self, other: "BasicGeometry") -> int:
        # return cmp(self.number, other.number)
        return (self.number > other.number) - (
            self.number < other.number
        )  # replacement for cmp function not availbale with Python3

    # ----------GENERATE PHASES------------------------------------------------

    def generate(self, rve: Rve = None) -> cq.Shape:
        """
        BasicGeometry class to manage shapes

        :param rve: Rve is required when the shape is a TPMS
        :type rve: Rve
        :return: CQ Shape
        """

        if self.shape.lower() == "box":
            cqshape = self.geometry.createBox()
        elif self.shape.lower() == "cylinder":
            cqshape = self.geometry.createCylinder()
        elif self.shape.lower() == "extrudedpolygon":
            cqshape = self.geometry.createExtrudedpolygon()
        elif self.shape.lower() == "capsule":
            cqshape = self.geometry.createCapsule()
        elif self.shape.lower() == "sphere":
            cqshape = self.geometry.createSphere()
        elif self.shape.lower() == "ellipsoid":
            cqshape = self.geometry.createEllipsoid()
        elif self.shape.lower() == "tpms":
            cqshape = self.geometry.createTpms(path_data=self.path_data, rve=rve)
        elif self.shape.lower() == "polyhedron":
            cqshape = self.geometry.createPolyhedron()
        else:
            raise ValueError(self.shape + " is not recognised")

        return cq.Shape(cqshape.val().wrapped)
