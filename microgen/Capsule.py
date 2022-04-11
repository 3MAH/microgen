from microgen.Functions import rotateEuler
import cadquery as cq

# ----------CAPSULE-----------------------------------------------------------------------------------------#


class Capsule:
    def __init__(self, center, angle, height, radius, number):
        self.center = center
        self.angle = angle
        self.height = height
        self.radius = radius
        self.number = number
        self.name_part = "capsule" + str(self.number)

    def createCapsule(self):
        cylinder = cq.Solid.makeCylinder(
            self.radius,
            self.height,
            pnt=cq.Vector(
                -self.height / 2.0 + self.center[0], self.center[1], self.center[2]
            ),
            dir=cq.Vector(0.1, 0.0, 0.0),
            angleDegrees=360,
        )
        sphereG = cq.Solid.makeSphere(
            self.radius,
            cq.Vector(
                self.center[0] - self.height / 2.0, self.center[1], self.center[2]
            ),
            angleDegrees1=-90,
        )
        sphereD = cq.Solid.makeSphere(
            self.radius,
            cq.Vector(
                self.center[0] + self.height / 2.0, self.center[1], self.center[2]
            ),
            angleDegrees1=-90,
        )
        capsule = cylinder.fuse(sphereG)
        capsule = capsule.fuse(sphereD)
        capsule = rotateEuler(capsule, self.center, self.angle[0], self.angle[1], self.angle[2])
        return cq.Workplane().add(capsule.Solids())
