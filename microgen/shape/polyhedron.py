"""Polyhedron.

=============================================
Polyhedron (:mod:`microgen.shape.polyhedron`)
=============================================
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Tuple

import cadquery as cq
import numpy as np
import pyvista as pv

from .shape import Shape

if TYPE_CHECKING:
    from microgen.shape import KwargsGenerateType, Vector3DType

Vertex = Tuple[float, float, float]
Face = Dict[str, List[int]]


class Polyhedron(Shape):
    """Class to generate a Polyhedron with a given set of faces and vertices.

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Polyhedron().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: Polyhedron,
        dic: dict[str, list[Vertex | Face]] | None = None,
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the polyhedron.

        .. warning:
            Give a center parameter only if the polyhedron must be translated
            from its original position.
        """
        super().__init__(**kwargs)
        if dic is None:
            self.dic: dict[str, list[Vertex | Face]] = {
                "vertices": [
                    (1.0, 1.0, 1.0),
                    (1.0, -1.0, -1.0),
                    (-1.0, 1.0, -1.0),
                    (-1.0, -1.0, 1.0),
                ],
                "faces": [
                    {"vertices": [0, 1, 2]},
                    {"vertices": [0, 3, 1]},
                    {"vertices": [0, 2, 3]},
                    {"vertices": [1, 2, 3]},
                ],
            }
        else:
            self.dic = dic

        self.faces_ixs = [face["vertices"] for face in self.dic["faces"]]
        for ixs in self.faces_ixs:
            ixs.append(ixs[0])

    def generate(self: Polyhedron, **_: KwargsGenerateType) -> cq.Shape:
        """Generate a polyhedron CAD shape using the given parameters."""
        faces = []
        for ixs in self.faces_ixs:
            lines = []
            for v1, v2 in zip(ixs, ixs[1:]):
                # tuple(map(sum, zip(a, b))) -> sum of tuples value by value
                vertice_coords1 = tuple(
                    map(sum, zip(self.center, self.dic["vertices"][v1])),
                )
                vertice_coords2 = tuple(
                    map(sum, zip(self.center, self.dic["vertices"][v2])),
                )
                lines.append(
                    cq.Edge.makeLine(
                        cq.Vector(*vertice_coords1),
                        cq.Vector(*vertice_coords2),
                    ),
                )
            wire = cq.Wire.assembleEdges(lines)
            faces.append(cq.Face.makeFromWires(wire))
        shell = cq.Shell.makeShell(faces)
        solid = cq.Solid.makeSolid(shell)
        return cq.Shape(solid.wrapped)

    def generate_vtk(self: Polyhedron, **_: KwargsGenerateType) -> pv.PolyData:
        """Generate a polyhedron VTK shape using the given parameters."""
        faces_pv = copy.deepcopy(self.faces_ixs)
        for vertices_in_face in faces_pv:
            del vertices_in_face[-1]
            vertices_in_face.insert(0, len(vertices_in_face))

        vertices = np.array(self.dic["vertices"])
        faces = np.hstack(faces_pv)

        return pv.PolyData(vertices, faces).compute_normals()

    def generateVtk(self: Polyhedron, **_: KwargsGenerateType) -> pv.PolyData:  # noqa: N802
        """Deprecated method. Use generate_vtk instead."""  # noqa: D401
        return self.generate_vtk()


def read_obj(filename: str) -> dict[str, list[Vertex | Face]]:
    """Read vertices and faces from obj format file for polyhedron."""
    dic: dict[str, list[Vertex | Face]] = {"vertices": [], "faces": []}
    with Path(filename).open(encoding="utf-8") as f:
        for line in f:
            data = line.split(" ")
            if data[0] == "v":
                dic["vertices"].append(tuple(float(x) for x in data[1:4]))
            elif data[0] == "f":
                dic["faces"].append(
                    {"vertices": [int(vertex) - 1 for vertex in data[1:]]},
                )
    return dic
