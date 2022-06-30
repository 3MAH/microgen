"""
=============================================
Polyhedron (:mod:`microgen.shape.polyhedron`)
=============================================
"""
import cadquery as cq
import pyvista as pv
import numpy as np
import copy

from .basicGeometry import BasicGeometry


class Polyhedron(BasicGeometry):
    """
    Class to generate a Polyhedron with a given set of faces and vertices
    """

    def __init__(
        self,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        dic: dict[str, list] = {
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
        },
    ) -> None:
        """
        .. warning:
            Give a center parameter only if the polyhedron must be translated from its original position.
        """
        super().__init__(shape="Polyhedron", center=center, orientation=orientation)
        self.dic = dic
        self.faces_ixs = [face["vertices"] for face in dic["faces"]]
        for ixs in self.faces_ixs:
            ixs.append(ixs[0])

    def generate(self) -> cq.Shape:
        faces = []
        for ixs in self.faces_ixs:
            lines = []
            for v1, v2 in zip(ixs, ixs[1:]):
                # tuple(map(sum, zip(a, b))) -> sum of tuples value by value
                vertice_coords1 = tuple(
                    map(sum, zip(self.center, self.dic["vertices"][v1]))
                )
                vertice_coords2 = tuple(
                    map(sum, zip(self.center, self.dic["vertices"][v2]))
                )
                lines.append(
                    cq.Edge.makeLine(
                        cq.Vector(*vertice_coords1), cq.Vector(*vertice_coords2)
                    )
                )
            wire = cq.Wire.assembleEdges(lines)
            faces.append(cq.Face.makeFromWires(wire))
        shell = cq.Shell.makeShell(faces)
        solid = cq.Solid.makeSolid(shell)
        return cq.Shape(solid.wrapped)

    def generateVtk(self) -> pv.PolyData:
        facesPv = copy.deepcopy(self.faces_ixs)
        for vertices_in_face in facesPv:
            del vertices_in_face[-1]
            vertices_in_face.insert(0, len(vertices_in_face))

        vertices = np.array(self.dic["vertices"])
        faces = np.hstack(facesPv)

        return pv.PolyData(vertices, faces)


def read_obj(filename: str):
    """
    Reads vertices and faces from obj format file for polyhedron
    """
    dic = {"vertices": [], "faces": []}
    with open(filename, "r") as f:
        for line in f:
            data = line.split(" ")
            if data[0] == "v":
                x = float(data[1])
                y = float(data[2])
                z = float(data[3])
                dic["vertices"].append([x, y, z])
            elif data[0] == "f":
                vertices = data[1:]
                for i in range(len(vertices)):
                    vertices[i] = int(vertices[i]) - 1
                dic["faces"].append({"vertices": vertices})
    return dic
