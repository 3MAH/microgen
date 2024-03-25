"""
=============================================
Sphere (:mod:`microgen.shape.sphere`)
=============================================
"""

from typing import Tuple

import cadquery as cq
import numpy as np
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

    def generate(self, **kwargs) -> cq.Shape:
        # Temporary workaround bug fix for OpenCascade bug using a random
        # direct parameter for cq.Workplane().sphere() method
        # Related to issue https://github.com/CadQuery/cadquery/issues/1461
        _seed = 38
        _random_direction_creation_axis = tuple(np.random.default_rng(_seed).random(3))
        sphere = (
            cq.Workplane()
            .sphere(radius=self.radius, direct=_random_direction_creation_axis)
            .translate(self.center)
        )
        return cq.Shape(sphere.val().wrapped)

    def generateVtk(
        self, theta_resolution=50, phi_resolution=50, **kwargs
    ) -> pv.PolyData:
        return pv.Sphere(
            radius=self.radius,
            center=tuple(self.center),
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        )
