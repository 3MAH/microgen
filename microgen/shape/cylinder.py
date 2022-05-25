from ..operations import rotateEuler
import cadquery as cq
import numpy as np

# ----------CYLINDER-----------------------------------------------------------------------------------------#


class Cylinder:
    def __init__(self, center: np.ndarray[float, float, float], angle: np.ndarray[float, float, float], height: float, radius: float, number: int) -> None:
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
