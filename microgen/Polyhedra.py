import cadquery as cq

# ----------POLYHEDRA-----------------------------------------------------------------------------------------#


class Polyhedra:
    def __init__(self, dic, n):
        self.dic = dic
        self.faces_ixs = [face["vertices"] for face in dic["faces"]]
        for ixs in self.faces_ixs:
            ixs.append(ixs[0])
        self.number = n
        self.name_part = "polyhedra" + str(self.number)

    def createPolyhedra(self):

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
