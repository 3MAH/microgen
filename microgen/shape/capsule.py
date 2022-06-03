import cadquery as cq
import numpy as np

import vtk
import pyvista as pv

from ..operations import rotateEuler
from ..pvoperations import rotatePvEuler

# ----------CAPSULE-------------------------------------------------------------------#


class Capsule:
    """
    Class to generate a capsule (cylinder with hemispherical ends)
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
            capsule, self.center, self.angle[0], self.angle[1], self.angle[2]
        )
        return cq.Workplane().add(capsule.Solids())

    def createPvCapsule(self, resolution=100, capping=True) -> pv.PolyData:
        capsule = vtk.vtkCapsuleSource()
        capsule.SetCenter(self.center[0], self.center[1], self.center[2])
        capsule.SetRadius(self.radius)
        capsule.SetCylinderLength(self.height)
        capsule.SetThetaResolution(resolution)
        capsule.SetPhiResolution(resolution)
        capsule.Update()
        capsule = pv.wrap(capsule.GetOutput())
        capsule = rotatePvEuler(capsule, self.center, self.angle[0], self.angle[1], self.angle[2])
        return capsule
