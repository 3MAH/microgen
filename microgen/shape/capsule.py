"""
=======================================
Capsule (:mod:`microgen.shape.capsule`)
=======================================
"""
import cadquery as cq

import pyvista as pv

from ..operations import rotateEuler, rotatePvEuler
from .basicGeometry import BasicGeometry


class Capsule(BasicGeometry):
    """
    Class to generate a capsule (cylinder with hemispherical ends)
    """

    def __init__(
        self,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        height: float = 1,
        radius: float = 0.5,
    ) -> None:
        super().__init__(shape="Capsule", center=center, orientation=orientation)
        self.height = height
        self.radius = radius

    def generate(self) -> cq.Shape:
        cylinder = cq.Solid.makeCylinder(
            self.radius,
            self.height,
            pnt=cq.Vector(
                -self.height / 2.0 + self.center[0], self.center[1], self.center[2]
            ),
            dir=cq.Vector(0.1, 0.0, 0.0),
            angleDegrees=360,
        )
        sphereL = cq.Solid.makeSphere(
            self.radius,
            cq.Vector(
                self.center[0] - self.height / 2.0, self.center[1], self.center[2]
            ),
            angleDegrees1=-90,
        )
        sphereR = cq.Solid.makeSphere(
            self.radius,
            cq.Vector(
                self.center[0] + self.height / 2.0, self.center[1], self.center[2]
            ),
            angleDegrees1=-90,
        )
        capsule = cylinder.fuse(sphereL)
        capsule = capsule.fuse(sphereR)
        capsule = rotateEuler(
            capsule,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return capsule

    def generateVtk(
        self, resolution=100, theta_resolution=50, phi_resolution=50, capping=True
    ) -> pv.PolyData:
        cylinder = pv.Cylinder(
            center=self.center,
            radius=self.radius,
            height=self.height,
            resolution=resolution,
            capping=capping,
        )
        sphereL = pv.Sphere(
            radius=self.radius,
            center=(self.center[0] - self.height / 2, self.center[1], self.center[2]),
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        )
        sphereR = pv.Sphere(
            radius=self.radius,
            center=(self.center[0] + self.height / 2, self.center[1], self.center[2]),
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        )
        capsule = cylinder.merge([sphereL, sphereR])
        capsule = rotatePvEuler(
            capsule,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return capsule
