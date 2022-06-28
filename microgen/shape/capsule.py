"""
=======================================
Capsule (:mod:`microgen.shape.capsule`)
=======================================
"""
import cadquery as cq

# import vtk
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
            capsule,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return capsule

    # def generateVtk(self, resolution=100, capping=True) -> pv.PolyData:
    #     capsule = vtk.vtkCapsuleSource()
    #     capsule.SetCenter(self.center[0], self.center[1], self.center[2])
    #     capsule.SetRadius(self.radius)
    #     capsule.SetCylinderLength(self.height)
    #     capsule.SetThetaResolution(resolution)
    #     capsule.SetPhiResolution(resolution)
    #     capsule.Update()
    #     capsule = pv.wrap(capsule.GetOutput())
    #     capsule = rotatePvEuler(capsule, self.center, 
    #                             self.orientation[0], self.orientation[1], self.orientation[2])
    #     return capsule
