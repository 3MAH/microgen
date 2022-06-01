import cadquery as cq
import pyvista as pv
import numpy as np

from ..operations import rotateEuler
from ..pvoperations import rotatePvEuler

# ----------CYLINDER-----------------------------------------------------------------------------------------#


class Cylinder:
    """
    Class to generate a cylinder
    """
    def __init__(
        self,
        center: np.ndarray,
        angle: np.ndarray,
        height: float,
        radius: float,
        number: int,
    ) -> None:
        self.center = center
        self.angle = angle
        self.radius = radius
        self.height = height
        self.number = number
        self.name_part = "cylinder" + str(self.number)

    def createCylinder(self) -> cq.Workplane:
        cylinder = (
            cq.Workplane("YZ")
            .circle(self.radius)
            .extrude(self.height)
            .translate(
                (self.center[0] - self.height / 2.0, self.center[1], self.center[2])
            )
        )
        cylinder = rotateEuler(
            cylinder, self.center, self.angle[0], self.angle[1], self.angle[2]
        )
        return cylinder


    def createPvCylinder(self, resolution=100, capping=True) -> pv.PolyData:
        cylinder = pv.Cylinder(
            center=tuple(self.center),
            direction=(1.0, 0.0, 0.0),
            radius=self.radius,
            height=self.height,
            resolution=resolution,
            capping=capping
        )
        cylinder = rotatePvEuler(
            cylinder, self.center, self.angle[0], self.angle[1], self.angle[2]
        )
        return cylinder
