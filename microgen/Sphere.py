from microgen.Functions import *
import numpy as np
import cadquery as cq

#----------SPHERE-----------------------------------------------------------------------------------------#
                
class sphere :
    def __init__(self,center,r,n):
        self.center=center
        self.radius=r
        self.number=n
        self.name_part='sphere' + str(self.number)
                
    def create_sphere(self):
        return cq.Workplane().sphere(self.radius).translate((self.center[0],self.center[1],self.center[2]))
