import cadquery as cq
import pyvista as pv
import numpy as np

from ..operations import rotateEuler
from ..pvoperations import rotatePvEuler

# ----------ELLIPSOID-----------------------------------------------------------------------------------------#


class Ellipsoid:
    """
    Class to generate an ellipsoid
    """
    def __init__(
        self,
        center: np.ndarray,
        angle: np.ndarray,
        a_x: float,
        a_y: float,
        a_z: float,
        number: int,
    ) -> None:
        self.center = center
        self.angle = angle
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
            ellipsoid, self.center, self.angle[0], self.angle[1], self.angle[2]
        )
        return cq.Workplane().add(ellipsoid)
        
    def createPvEllipsoid(self) -> pv.PolyData:
        transform_matrix = np.array(
            [
                [self.a_x, 0, 0, self.center[0]],
                [0, self.a_y, 0, self.center[1]],
                [0, 0, self.a_z, self.center[2]],
            ]
        )

        sphere = pv.Sphere(
            radius=self.radius,
            center=tuple(self.center)
        )
        ellipsoid = sphere.transform(transform_matrix)
        ellipsoid = rotatePvEuler(
            ellipsoid, self.center, self.angle[0], self.angle[1], self.angle[2]
        )
        return ellipsoid
