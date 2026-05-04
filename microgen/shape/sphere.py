"""Sphere.

=====================================
Sphere (:mod:`microgen.shape.sphere`)
=====================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pyvista as pv

from .shape import Shape

if TYPE_CHECKING:
    from microgen.cad import CadShape
    from microgen.shape import KwargsGenerateType, Vector3DType


class Sphere(Shape):
    """Class to generate a sphere.

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Sphere().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: Sphere,
        radius: float = 1,
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the sphere."""
        super().__init__(**kwargs)
        self.radius = radius

    def generate(self: Sphere, **_: KwargsGenerateType) -> CadShape:
        """Generate a sphere CAD shape (OCCT).  Requires the ``[cad]`` extra."""
        from microgen.cad import make_sphere

        return make_sphere(self.radius, self.center)

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
