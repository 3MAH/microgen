"""Extruded Polygon.

========================================================
Extruded Polygon (:mod:`microgen.shape.extruded_polygon`)
========================================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Sequence

import cadquery as cq
import numpy as np
import pyvista as pv

from microgen.operations import rotateEuler, rotatePvEuler

from .shape import Shape

if TYPE_CHECKING:
    from microgen.shape import KwargsGenerateType, Vector3DType


class ExtrudedPolygon(Shape):
    """ExtrudedPolygon.

    Class to generate an extruded polygon with a given list of points and a thickness

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.ExtrudedPolygon().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: ExtrudedPolygon,
        list_corners: Sequence[tuple[float, float]] | None = None,
        height: float = 1,
        **kwargs: Vector3DType | dict[str, Sequence[tuple[float, float]]],
    ) -> None:
        """Initialize the extruded polygon. The default shape is a hexagon."""
        center = kwargs.get("center", (0, 0, 0))
        orientation = kwargs.get("orientation", (0, 0, 0))
        super().__init__(center=center, orientation=orientation)

        if kwargs.get("listCorners", None) is not None:
            list_corners = kwargs["listCorners"]

        if list_corners is None:
            self.list_corners = [
                (1, 0),
                (0.5, 0.5 * np.sqrt(3)),
                (-0.5, 0.5 * np.sqrt(3)),
                (-1, 0),
                (-0.5, -0.5 * np.sqrt(3)),
                (0.5, -0.5 * np.sqrt(3)),
                (1, 0),
            ]  # hexagon
        else:
            self.list_corners = list_corners
        self.height = height

    def generate(self: ExtrudedPolygon, **_: KwargsGenerateType) -> cq.Shape:
        """Generate an extruded polygon CAD shape using the given parameters."""
        poly = (
            cq.Workplane("YZ")
            .polyline(self.list_corners)
            .close()
            .extrude(self.height)
            .translate(
                (self.center[0] - self.height / 2.0, self.center[1], self.center[2]),
            )
        )
        poly = rotateEuler(
            poly,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return cq.Shape(poly.val().wrapped)

    def generate_vtk(self: ExtrudedPolygon, **_: KwargsGenerateType) -> pv.PolyData:
        """Generate an extruded polygon VTK shape using the given parameters."""
        vertices = [
            [
                self.center[0] + 0.5 * self.height,
                self.center[1] + corner[0],
                self.center[2] + corner[1],
            ]
            for corner in self.list_corners
        ]
        faces = np.arange(len(vertices))
        faces = np.insert(faces, 0, len(vertices))

        poly = (
            pv.PolyData(vertices, faces)
            .extrude([self.height, 0, 0], capping=True)
            .compute_normals()
        )

        return rotatePvEuler(
            poly,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )

    def generateVtk(self: ExtrudedPolygon, **_: KwargsGenerateType) -> pv.PolyData:  # noqa: N802
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk()
