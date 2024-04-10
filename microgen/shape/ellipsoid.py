"""Ellipsoid.

=============================================
Ellipsoid (:mod:`microgen.shape.ellipsoid`)
=============================================
"""

from __future__ import annotations

import cadquery as cq
import numpy as np
import pyvista as pv

from microgen.operations import rotateEuler, rotatePvEuler

from .basic_geometry import BasicGeometry


class Ellipsoid(BasicGeometry):
    """Class to generate an ellipsoid.

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Ellipsoid().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        radii: tuple[float, float, float] = (1, 0.5, 0.25),
    ) -> None:
        """Initialize the ellipsoid."""
        super().__init__(shape="Ellipsoid", center=center, orientation=orientation)
        self.radii = radii

    def generate(self, **_) -> cq.Shape:
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

    def generate_vtk(self, **_) -> pv.PolyData:
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

    def generateVtk(self, **kwargs) -> pv.PolyData:  # noqa: N802
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk(**kwargs)
