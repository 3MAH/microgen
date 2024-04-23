"""Cylinder.

=========================================
Cylinder (:mod:`microgen.shape.cylinder`)
=========================================
"""

from __future__ import annotations

from typing import Any

import cadquery as cq
import pyvista as pv

from microgen.operations import rotateEuler, rotatePvEuler

from .basic_geometry import BasicGeometry


class Cylinder(BasicGeometry):
    """Class to generate a cylinder.

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Cylinder().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: Cylinder,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        height: float = 1,
        radius: float = 0.5,
    ) -> None:
        """Initialize the cylinder."""
        super().__init__(shape="Cylinder", center=center, orientation=orientation)
        self.radius = radius
        self.height = height

    def generate(self: Cylinder, **_: dict[str, Any]) -> cq.Shape:
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
        **_: dict[str, Any],
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
        **kwargs: dict[str, Any],
    ) -> pv.PolyData:
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk(
            resolution=resolution,
            **kwargs,
        )
