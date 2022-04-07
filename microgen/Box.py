from microgen.Functions import rotateEuler
import cadquery as cq

# ----------BOX-----------------------------------------------------------------------------------------#
# MB 03/12/2021


class Box:
    def __init__(self, center, angle, a1, a2, a3, number):
        self.center = center
        self.angle = angle
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.number = number
        self.name_part = "box" + str(self.number)

    def createBox(self):
        box = (
            cq.Workplane()
            .box(self.a1, self.a2, self.a3)
            .translate((self.center[0], self.center[1], self.center[2]))
        )
        box = rotateEuler(box, self.center, self.angle[0], self.angle[1], self.angle[2])
        return box
