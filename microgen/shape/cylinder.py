"""
Cylinder.

=========================================
Cylinder (:mod:`microgen.shape.cylinder`)
=========================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pyvista as pv

from microgen.operations import rotate

from .shape import Shape

if TYPE_CHECKING:
    from microgen.cad import CadShape
    from microgen.shape import KwargsGenerateType, Vector3DType


class Cylinder(Shape):
    """
    Class to generate a cylinder.

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Cylinder().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: Cylinder,
        height: float = 1,
        radius: float = 0.5,
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the cylinder."""
        super().__init__(**kwargs)
        self.radius = radius
        self.height = height

    def generate(self: Cylinder, **_: KwargsGenerateType) -> CadShape:
        """Generate a cylinder CAD shape (OCCT).  Requires the ``[cad]`` extra."""
        from microgen.cad import make_cylinder  # noqa: PLC0415

        shape = make_cylinder(
            radius=self.radius,
            height=self.height,
            center=self.center,
            axis=(1.0, 0.0, 0.0),
        )
        return rotate(shape, self.center, self.orientation)

    def generate_vtk(
        self: Cylinder,
        resolution: int = 100,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate a cylinder VTK shape using the given parameters."""
        cylinder = pv.Cylinder(
            center=tuple(self.center),
            direction=(1.0, 0.0, 0.0),
            radius=self.radius,
            height=self.height,
            resolution=resolution,
            capping=True,
        )
        return rotate(cylinder, self.center, self.orientation)

    def generateVtk(  # noqa: N802
        self: Cylinder,
        resolution: int = 100,
        **kwargs: KwargsGenerateType,
    ) -> pv.PolyData:
        """Deprecated. Use :meth:`generate_vtk` instead."""
        return self.generate_vtk(
            resolution=resolution,
            **kwargs,
        )
