"""Capsule.

=======================================
Capsule (:mod:`microgen.shape.capsule`)
=======================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import cadquery as cq
import pyvista as pv

from microgen.operations import rotateEuler, rotatePvEuler

from .shape import Shape

if TYPE_CHECKING:
    from microgen.shape import KwargsGenerateType, Vector3DType


class Capsule(Shape):
    """Class to generate a capsule (cylinder with hemispherical ends).

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Capsule().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: Capsule,
        height: float = 1,
        radius: float = 0.5,
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the capsule."""
        super().__init__(**kwargs)
        self.height = height
        self.radius = radius

    def generate(self: Capsule, **_: KwargsGenerateType) -> cq.Shape:
        """Generate a capsule CAD shape using the given parameters."""
        cylinder = cq.Solid.makeCylinder(
            self.radius,
            self.height,
            pnt=cq.Vector(
                -self.height / 2.0 + self.center[0],
                self.center[1],
                self.center[2],
            ),
            dir=cq.Vector(0.1, 0.0, 0.0),
            angleDegrees=360,
        )
        sphere_left = cq.Solid.makeSphere(
            self.radius,
            cq.Vector(
                self.center[0] - self.height / 2.0,
                self.center[1],
                self.center[2],
            ),
            angleDegrees1=-90,
        )
        sphere_right = cq.Solid.makeSphere(
            self.radius,
            cq.Vector(
                self.center[0] + self.height / 2.0,
                self.center[1],
                self.center[2],
            ),
            angleDegrees1=-90,
        )
        capsule = cylinder.fuse(sphere_left)
        capsule = capsule.fuse(sphere_right)
        return rotateEuler(
            capsule,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )

    def generate_vtk(
        self: Capsule,
        resolution: int = 100,
        theta_resolution: int = 50,
        phi_resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate a capsule VTK shape using the given parameters."""
        cylinder = pv.Cylinder(
            center=self.center,
            radius=self.radius,
            height=self.height,
            resolution=resolution,
            capping=True,
        ).triangulate()
        sphere_left = pv.Sphere(
            radius=self.radius,
            center=(self.center[0] - self.height / 2, self.center[1], self.center[2]),
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        ).triangulate()
        sphere_right = pv.Sphere(
            radius=self.radius,
            center=(self.center[0] + self.height / 2, self.center[1], self.center[2]),
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        ).triangulate()
        capsule = cylinder.boolean_union(sphere_left).boolean_union(sphere_right)
        return rotatePvEuler(
            capsule,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )

    def generateVtk(  # noqa: N802
        self: Capsule,
        resolution: int = 100,
        theta_resolution: int = 50,
        phi_resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk(
            resolution=resolution,
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        )
