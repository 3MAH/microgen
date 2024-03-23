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
        radii: Tuple[float, float, float] = (1, 0.5, 0.25),
    ) -> None:
        super().__init__(shape="Ellipsoid", center=center, orientation=orientation)
        self.radii = radii

    def generate(self) -> cq.Shape:
        transform_mat = cq.Matrix(
            [
                [self.radii[0], 0, 0, self.center[0]],
                [0, self.radii[1], 0, self.center[1]],
                [0, 0, self.radii[2], self.center[2]],
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

    def generateVtk(self) -> pv.PolyData:
        transform_matrix = np.array(
            [
                [self.radii[0], 0, 0, self.center[0]],
                [0, self.radii[1], 0, self.center[1]],
                [0, 0, self.radii[2], self.center[2]],
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
