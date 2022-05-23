import cadquery as cq

# ----------SPHERE--------------------------------------------------------#


class Sphere:
    def __init__(self, center, radius, number):
        self.center = center
        self.radius = radius
        self.number = number
        self.name_part = "sphere" + str(self.number)

    def createSphere(self):
        return (
            cq.Workplane()
            .sphere(self.radius)
            .translate((self.center[0], self.center[1], self.center[2]))
        )
