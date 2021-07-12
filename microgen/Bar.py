from microgen.Functions import *
import numpy as np
import cadquery as cq

#----------BAR-----------------------------------------------------------------------------------------#
                
class bar :
    def __init__(self,center,angle,h,r,n):
        self.center=center
        self.angle=angle
        self.radius=r
        self.height=h
        self.number=n
        self.name_part='bar' + str(self.number)
                
    def create_bar(self):
        cylinder = cq.Solid.makeCylinder(self.radius, self.height, pnt = cq.Vector(-self.height/2.0+self.center[0], self.center[1], self.center[2]), dir = cq.Vector(0.1, 0.0, 0.0), angleDegrees=360)
        sphereG = cq.Solid.makeSphere(self.radius, cq.Vector(self.center[0]-self.height/2.0,self.center[1],self.center[2]), angleDegrees1=-90)
        sphereD = cq.Solid.makeSphere(self.radius, cq.Vector(self.center[0]+self.height/2.0,self.center[1],self.center[2]), angleDegrees1=-90)
        bar = cq.Compound.makeCompound([cylinder_bar,])
        bar = bar.fuse(sphereG_bar)
        bar = bar.fuse(sphereD_bar)
        bar = rotateEuler(bar, self.center, self.angle[0], self.angle[1], self.angle[2])
        return bar
		
