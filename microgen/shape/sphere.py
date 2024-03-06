"""
=============================================
Sphere (:mod:`microgen.shape.sphere`)
=============================================
"""
from typing import Tuple

import cadquery as cq
import pyvista as pv

from .basicGeometry import BasicGeometry


class Sphere(BasicGeometry):
    """
    Class to generate a sphere

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Sphere().generateVtk()
       shape.plot(color='white')
    """

    def __init__(
        self,
        center: Tuple[float, float, float] = (0, 0, 0),
        radius: float = 1,
    ) -> None:
        super().__init__(shape="Sphere", center=center)
        self.radius = radius

    def generate(self) -> cq.Shape:
        sphere = (
            cq.Workplane()
            .sphere(self.radius)
            .translate((self.center[0], self.center[1], self.center[2]))
        )
        return cq.Shape(sphere.val().wrapped)

    def generateVtk(self, theta_resolution=30, phi_resolution=30) -> pv.PolyData:
        return pv.Sphere(
            radius=self.radius,
            center=tuple(self.center),
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        )
