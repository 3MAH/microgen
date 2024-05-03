"""Cylinder.

=========================================
Cylinder (:mod:`microgen.shape.cylinder`)
=========================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import cadquery as cq
import pyvista as pv

from microgen.operations import rotateEuler, rotatePvEuler

from .shape import Shape

if TYPE_CHECKING:
    from microgen.shape import KwargsGenerateType, Vector3DType


class Cylinder(Shape):
    """Class to generate a cylinder.

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

    def generate(self: Cylinder, **_: KwargsGenerateType) -> cq.Shape:
        """Generate a cylinder CAD shape using the given parameters."""
        cylinder = (
            cq.Workplane("YZ")
            .circle(self.radius)
            .extrude(self.height)
            .translate(
                (self.center[0] - self.height / 2.0, self.center[1], self.center[2]),
            )
        )
        cylinder = rotateEuler(
            cylinder,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return cq.Shape(cylinder.val().wrapped)

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
        return rotatePvEuler(
            cylinder,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )

    def generateVtk(  # noqa: N802
        self: Cylinder,
        resolution: int = 100,
        **kwargs: KwargsGenerateType,
    ) -> pv.PolyData:
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk(
            resolution=resolution,
            **kwargs,
        )
