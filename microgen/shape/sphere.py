import cadquery as cq
import pyvista as pv
import numpy as np

# ----------SPHERE--------------------------------------------------------#


class Sphere:
    """
    Class to generate a sphere
    """
    def __init__(self, center: np.ndarray, radius: float, number: int) -> None:
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

    def createPvSphere(self, theta_resolution=30, phi_resolution=30) -> pv.PolyData:
        return pv.Sphere(
            radius=self.radius,
            center=tuple(self.center)
        )
