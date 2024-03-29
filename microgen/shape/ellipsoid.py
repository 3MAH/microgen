"""
=============================================
Ellipsoid (:mod:`microgen.shape.ellipsoid`)
=============================================
"""

from typing import Tuple

import cadquery as cq
import numpy as np
import pyvista as pv

from ..operations import rotateEuler, rotatePvEuler
from .basicGeometry import BasicGeometry


class Ellipsoid(BasicGeometry):
    """
    Class to generate an ellipsoid

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Ellipsoid().generateVtk()
       shape.plot(color='white')
    """

    def __init__(
        self,
        center: Tuple[float, float, float] = (0, 0, 0),
        orientation: Tuple[float, float, float] = (0, 0, 0),
        a_x: float = 1,
        a_y: float = 0.5,
        a_z: float = 0.25,
    ) -> None:
        super().__init__(shape="Ellipsoid", center=center, orientation=orientation)
        self.a_x = a_x
        self.a_y = a_y
        self.a_z = a_z

    def generate(self, **kwargs) -> cq.Shape:
        transform_mat = cq.Matrix(
            [
                [self.a_x, 0, 0, self.center[0]],
                [0, self.a_y, 0, self.center[1]],
                [0, 0, self.a_z, self.center[2]],
            ]
        )

        sphere = cq.Solid.makeSphere(1.0, cq.Vector(0, 0, 0), angleDegrees1=-90)
        ellipsoid = sphere.transformGeometry(transform_mat)
        ellipsoid = rotateEuler(
            ellipsoid,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return ellipsoid

    def generateVtk(self, **kwargs) -> pv.PolyData:
        transform_matrix = np.array(
            [
                [self.a_x, 0, 0, self.center[0]],
                [0, self.a_y, 0, self.center[1]],
                [0, 0, self.a_z, self.center[2]],
                [0, 0, 0, 1],
            ]
        )
        sphere = pv.Sphere(radius=1)
        ellipsoid = sphere.transform(transform_matrix, inplace=False)
        ellipsoid = rotatePvEuler(
            ellipsoid,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return ellipsoid
