"""
=============================================
Ellipsoid (:mod:`microgen.shape.ellipsoid`)
=============================================
"""
import cadquery as cq
import numpy as np

from ..operations import rotateEuler

# ----------ELLIPSOID-----------------------------------------------------------------------------------------#


class Ellipsoid:
    """
    Class to generate an ellipsoid
    """
    def __init__(
        self,
        center: np.ndarray,
        orientation: np.ndarray,
        a_x: float,
        a_y: float,
        a_z: float,
        number: int = 0,
    ) -> None:
        self.center = center
        self.orientation = orientation
        self.a_x = a_x
        self.a_y = a_y
        self.a_z = a_z
        self.number = number
        self.name_part = "ellipsoid" + str(self.number)

    def createEllipsoid(self) -> cq.Workplane:
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
            ellipsoid, self.center, self.orientation[0], self.orientation[1], self.orientation[2]
        )
        return cq.Workplane().add(ellipsoid)
