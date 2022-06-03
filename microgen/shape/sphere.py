"""
=============================================
Sphere (:mod:`microgen.shape.sphere`)
=============================================
"""
import cadquery as cq
import numpy as np


class Sphere:
    """
    Class to generate a sphere
    """
    def __init__(self, center: np.ndarray, radius: float, number: int = 0) -> None:
        self.center = center
        self.radius = radius
        self.number = number
        self.name_part = "sphere" + str(self.number)

    def createSphere(self) -> cq.Workplane:
        return (
            cq.Workplane()
            .sphere(self.radius)
            .translate((self.center[0], self.center[1], self.center[2]))
        )
