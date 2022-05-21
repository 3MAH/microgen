import numpy as np
import cadquery as cq

from microgen.Box import Box
from microgen.Sphere import Sphere
from microgen.Cylinder import Cylinder
from microgen.ExtrudedPolygon import ExtrudedPolygon
from microgen.Ellipsoid import Ellipsoid
from microgen.Capsule import Capsule
from microgen.Tpms import Tpms
from microgen.Polyhedron import Polyhedron


class BasicGeometry:
    def __init__(
        self,
        shape,
        param_geom,
        number=0,
        xc=0,
        yc=0,
        zc=0,
        psi=0,
        theta=0,
        phi=0,
        path_data=None,
    ):
        """DESCRIPTION

        Parameters
        ----------
        number : TYPE
            DESCRIPTION
        shape : TYPE
            DESCRIPTION
        xc, yc, zc : TYPE
            DESCRIPTION
        psi, theta, phi : TYPE
            DESCRIPTION
        param_geom : TYPE
            DESCRIPTION
        path_data : TYPE
            DESCRIPTION
        """
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
        self.angle = np.array([self.psi, self.theta, self.phi])
        self.name = self.shape + str(self.number)

        if self.shape.lower() == "box":
            self.geometry = Box(
                center=self.center,
                angle=self.angle,
                dim_x=self.param_geom["dim_x"],
                dim_y=self.param_geom["dim_y"],
                dim_z=self.param_geom["dim_z"],
                number=self.number,
            )
        if self.shape.lower() == "cylinder":
            self.geometry = Cylinder(
                center=self.center,
                angle=self.angle,
                height=self.param_geom["height"],
                radius=self.param_geom["radius"],
                number=self.number,
            )
        if self.shape.lower() == "extrudedpolygon":
            self.geometry = ExtrudedPolygon(
                center=self.center,
                angle=self.angle,
                listCorners=self.param_geom["listCorners"],
                height=self.param_geom["height"],
                number=self.number,
            )
        if self.shape.lower() == "capsule":
            self.geometry = Capsule(
                center=self.center,
                angle=self.angle,
                height=self.param_geom["height"],
                radius=self.param_geom["radius"],
                number=self.number,
            )
        if self.shape.lower() == "sphere":
            self.geometry = Sphere(
                center=self.center,
                radius=self.param_geom["radius"],
                number=self.number
            )
        if self.shape.lower() == "ellipsoid":
            self.geometry = Ellipsoid(
                center=self.center,
                angle=self.angle,
                a_x=self.param_geom["a_x"],
                a_y=self.param_geom["a_y"],
                a_z=self.param_geom["a_z"],
                number=self.number,
            )
        if self.shape.lower() == "tpms":
            if self.param_geom["type_surface"] == "custom":
                self.geometry = Tpms(
                    center=self.center,
                    angle=self.angle,
                    type_surface=self.param_geom["type_surface"],
                    type_part=self.param_geom["type_part"],
                    thickness=self.param_geom["thickness"],
                    number=self.number,
                    function=self.param_geom["function"],
                )
            else:
                self.geometry = Tpms(
                    center=self.center,
                    angle=self.angle,
                    type_surface=self.param_geom["type_surface"],
                    type_part=self.param_geom["type_part"],
                    thickness=self.param_geom["thickness"],
                    number=self.number,
                )

        if self.shape.lower() == "polyhedron":
            self.geometry = Polyhedron(dic=self.param_geom["dic"], 
                                       number=self.number)

    def __cmp__(self, other):
        # return cmp(self.number, other.number)
        return (self.number > other.number) - (
            self.number < other.number
        )  # replacement for cmp function not availbale with Python3

    # ----------GENERATE PHASES------------------------------------------------

    def generate(self, rve=None):
        """DESCRIPTION

        Parameters
        ----------
        rve : TYPE, optional
            DESCRIPTION

        Returns
        -------
        cq.Shape(cqshape.val().wrapped) : TYPE
            DESCRIPTION
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
            cqshape = self.geometry.createTpms(path_data=self.path_data,
                                               rve=rve)
        elif self.shape.lower() == "polyhedron":
            cqshape = self.geometry.createPolyhedron()
        else:
            raise ValueError(self.shape + " is not recognised")

        return cq.Shape(cqshape.val().wrapped)
