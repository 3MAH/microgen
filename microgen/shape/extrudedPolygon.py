from ..operations import rotateEuler
import cadquery as cq

# ----------ExtrudedPolygon-----------------------------------------------------------------------------------------#


class ExtrudedPolygon:
    def __init__(self, center, angle, listCorners, height, number):
        self.center = center
        self.angle = angle
        self.listCorners = listCorners
        self.height = height
        self.number = number
        self.name_part = "extrudedpolygon" + str(self.number)

    def createExtrudedpolygon(self):
        poly = (
            cq.Workplane("YZ")
            .polyline(self.listCorners)
            .close()
            .extrude(self.height)
            .translate(
                (self.center[0] - self.height / 2.0, self.center[1], self.center[2])
            )
        )
        poly = rotateEuler(
            poly, self.center, self.angle[0], self.angle[1], self.angle[2]
        )
        return poly
