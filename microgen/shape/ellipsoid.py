"""Ellipsoid.

=============================================
Ellipsoid (:mod:`microgen.shape.ellipsoid`)
=============================================
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import cadquery as cq
import numpy as np
import pyvista as pv

from microgen.operations import rotateEuler, rotatePvEuler

from .shape import Shape

if TYPE_CHECKING:
    from microgen.shape import KwargsGenerateType, Vector3DType


class Ellipsoid(Shape):
    """Class to generate an ellipsoid.

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Ellipsoid().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: Ellipsoid,
        radii: tuple[float, float, float] = (1, 0.5, 0.25),
        a_x: float | None = None,
        a_y: float | None = None,
        a_z: float | None = None,
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the ellipsoid."""
        super().__init__(**kwargs)
        if a_x is not None or a_y is not None or a_z is not None:
            warnings.warn(
                "The 'a_x', 'a_y', and 'a_z' parameters are deprecated. \
                    Use 'radii' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            if a_x is None:
                a_x = radii[0]
            if a_y is None:
                a_y = radii[1]
            if a_z is None:
                a_z = radii[2]
            radii = (a_x, a_y, a_z)

        self.radii = radii

    def generate(self: Ellipsoid, **_: KwargsGenerateType) -> cq.Shape:
        """Generate an ellipsoid CAD shape using the given parameters."""
        transform_mat = cq.Matrix(
            [
                [self.radii[0], 0, 0, self.center[0]],
                [0, self.radii[1], 0, self.center[1]],
                [0, 0, self.radii[2], self.center[2]],
            ],
        )

        sphere = cq.Solid.makeSphere(1.0, cq.Vector(0, 0, 0), angleDegrees1=-90)
        ellipsoid = sphere.transformGeometry(transform_mat)
        return rotateEuler(
            ellipsoid,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )

    def generate_vtk(self: Ellipsoid, **_: KwargsGenerateType) -> pv.PolyData:
        """Generate an ellipsoid VTK polydta using the given parameters."""
        transform_matrix = np.array(
            [
                [self.radii[0], 0, 0, self.center[0]],
                [0, self.radii[1], 0, self.center[1]],
                [0, 0, self.radii[2], self.center[2]],
                [0, 0, 0, 1],
            ],
        )
        sphere = pv.Sphere(radius=1)
        ellipsoid = sphere.transform(transform_matrix, inplace=False)
        return rotatePvEuler(
            ellipsoid,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )

    def generateVtk(self: Ellipsoid, **_: KwargsGenerateType) -> pv.PolyData:  # noqa: N802
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk()
