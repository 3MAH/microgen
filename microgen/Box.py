from microgen.Functions import *
import numpy as np
import cadquery as cq

#----------BOX-----------------------------------------------------------------------------------------#
# MB 03/12/2021
                
class box :
    def __init__(self,center,angle,a1,a2,a3,n):
        self.center=center
        self.angle=angle
        self.a1=a1
        self.a2=a2
        self.a3=a3
        self.number=n
        self.name_part='box' + str(self.number)
                
    def create_box(self):
        box = cq.Workplane().box(self.a1,self.a2,self.a3).translate((self.center[0],self.center[1],self.center[2]))
        box = rotateEuler(box, self.center, self.angle[0], self.angle[1], self.angle[2])
        return box


    
