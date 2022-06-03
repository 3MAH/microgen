import cadquery as cq
import pyvista as pv
import numpy as np
import copy

# ----------POLYHEDRON-----------------------------------------------------------------------------------------#


class Polyhedron:
    """
    Class to generate a Polyhedron with a given set of faces and vertices
    """
    def __init__(self, dic: dict[str, list], number: int) -> None:
        self.dic = dic
        self.faces_ixs = [face["vertices"] for face in dic["faces"]]
        for ixs in self.faces_ixs:
            ixs.append(ixs[0])
        self.number = number
        self.name_part = "polyhedron" + str(self.number)

    def createPolyhedron(self) -> cq.Workplane:

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
        return cq.Workplane().add(solid)

    def createPvPolyhedron(self) -> pv.PolyData:

        facesPv = copy.deepcopy(self.faces_ixs)
        for vertices_in_face in facesPv:
            del vertices_in_face[-1]
            vertices_in_face.insert(0,len(vertices_in_face))

        vertices = np.array(self.dic["vertices"])
        faces = np.hstack(facesPv)

        return pv.PolyData(vertices, faces)



        
