from microgen.Functions import *
import numpy as np
import cadquery as cq
#----------ELLIPSOID-----------------------------------------------------------------------------------------#
                
class ellipsoid :
    def __init__(self,center,angle,a1,a2,a3):
        self.center=center
        self.angle=angle
        self.a1=a1
        self.a2=a2
        self.a3=a3
        self.number=n
        self.name_part='ellipsoid' + str(self.number)
                
    def create_ellipsoid(self):
        transform_mat = cq.Matrix([[a1, 0, 0, self.center[0]],
                 [0, a2, 0, self.center[1]],
                 [0, 0, a3, self.center[2]]])
            
        sphere = cq.Solid.makeSphere(1.0, cq.Vector(0,0,0), angleDegrees1=-90)
        ellipsoid = sphere.transformGeometry(transform_mat)
        ellipsoid = rotateEuler(ellipsoid, self.center, self.angle[0], self.angle[1], self.angle[2])
        return ellipsoid
