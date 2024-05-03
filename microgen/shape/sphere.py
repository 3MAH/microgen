"""Sphere.

=====================================
Sphere (:mod:`microgen.shape.sphere`)
=====================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import cadquery as cq
import numpy as np
import pyvista as pv

from .basic_geometry import BasicGeometry

if TYPE_CHECKING:
    from microgen.shape import KwargsGenerateType, Vector3DType


class Sphere(BasicGeometry):
    """Class to generate a sphere.

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Sphere().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: Sphere,
        center: Vector3DType = (0, 0, 0),
        radius: float = 1,
    ) -> None:
        """Initialize the sphere."""
        super().__init__(shape="Sphere", center=center)
        self.radius = radius

    def generate(self: Sphere, **_: KwargsGenerateType) -> cq.Shape:
        """Generate a sphere CAD shape using the given parameters."""
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

    def generate_vtk(
        self: Sphere,
        theta_resolution: int = 50,
        phi_resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate a sphere VTK shape using the given parameters."""
        return pv.Sphere(
            radius=self.radius,
            center=tuple(self.center),
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        )

    def generateVtk(  # noqa: N802
        self: Sphere,
        theta_resolution: int = 50,
        phi_resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Deprecated method. Use generate_vtk instead."""  # noqa: D401
        return self.generate_vtk(theta_resolution, phi_resolution)
