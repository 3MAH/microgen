"""
=======================================
Capsule (:mod:`microgen.shape.capsule`)
=======================================
"""
import cadquery as cq
import numpy as np

from ..operations import rotateEuler


class Capsule:
    """
    Class to generate a capsule (cylinder with hemispherical ends)
    """
    def __init__(
        self,
        center: np.ndarray,
        orientation: np.ndarray,
        height: float,
        radius: float,
        number: int = 0,
    ) -> None:
        self.center = center
        self.orientation = orientation
        self.height = height
        self.radius = radius
        self.number = number
        self.name_part = "capsule" + str(self.number)

    def createCapsule(self) -> cq.Workplane:
        cylinder = cq.Solid.makeCylinder(
            self.radius,
            self.height,
            pnt=cq.Vector(
                -self.height / 2.0 + self.center[0], self.center[1], self.center[2]
            ),
            dir=cq.Vector(0.1, 0.0, 0.0),
            angleDegrees=360,
        )
        sphereG = cq.Solid.makeSphere(
            self.radius,
            cq.Vector(
                self.center[0] - self.height / 2.0, self.center[1], self.center[2]
            ),
            angleDegrees1=-90,
        )
        sphereD = cq.Solid.makeSphere(
            self.radius,
            cq.Vector(
                self.center[0] + self.height / 2.0, self.center[1], self.center[2]
            ),
            angleDegrees1=-90,
        )
        capsule = cylinder.fuse(sphereG)
        capsule = capsule.fuse(sphereD)
        capsule = rotateEuler(
            capsule, self.center, self.orientation[0], self.orientation[1], self.orientation[2]
        )
        return cq.Workplane().add(capsule.Solids())
