"""
=============================================
Polyhedron (:mod:`microgen.shape.polyhedron`)
=============================================
"""
import cadquery as cq

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
            "vertices": [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)],
            "faces": {"vertices": [(0, 1, 2), (0, 3, 1), (0, 2, 3), (1, 2, 3)]},
        },
    ) -> None:
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
                vertice_coords1 = self.dic["vertices"][v1]
                vertice_coords2 = self.dic["vertices"][v2]
                lines.append(
                    cq.Edge.makeLine(
                        cq.Vector(*vertice_coords1), cq.Vector(*vertice_coords2)
                    )
                )
            wire = cq.Wire.assembleEdges(lines)
            faces.append(cq.Face.makeFromWires(wire))
        shell = cq.Shell.makeShell(faces)
        solid = cq.Solid.makeSolid(shell)
        # return cq.Shape(cq.Workplane().add(solid).val().wrapped)  # need better way to convert from cq.Solid to cq.Shape
        return cq.Shape(solid.wrapped)


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
