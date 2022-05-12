from microgen.Operations import rotateEuler
import cadquery as cq

# ----------ELLIPSOID-----------------------------------------------------------------------------------------#


class Ellipsoid:
    def __init__(self, center, angle, a_x, a_y, a_z, number):
        self.center = center
        self.angle = angle
        self.a_x = a_x
        self.a_y = a_y
        self.a_z = a_z
        self.number = number
        self.name_part = "ellipsoid" + str(self.number)

    def createEllipsoid(self):
        transform_mat = cq.Matrix(
            [
                [self.a_x, 0, 0, self.center[0]],
                [0, self.a_y, 0, self.center[1]],
                [0, 0, self.a_z, self.center[2]],
            ]
        )

        sphere = cq.Solid.makeSphere(1.0, cq.Vector(0, 0, 0), angleDegrees1=-90)
        ellipsoid = sphere.transformGeometry(transform_mat)
        ellipsoid = rotateEuler(
            ellipsoid, self.center, self.angle[0], self.angle[1], self.angle[2]
        )
        return cq.Workplane().add(ellipsoid)
