import os
import sys
import numpy as np
import cadquery as cq
import gmsh
from microgen import *

#Size of the mesh
size_mesh = 0.03
a = 1.0
b = 1.0
c = 1.0

periodicity = 0
Revel = Rve(a,b,c,size_mesh)
phases = []

elem = BasicGeometry(0,'tpms',0.5,0.5,0.5,0.,0.,0.,['sheet','na'],'data')
elem.geometry.createSurfaces('gyroid', rve=Revel, thickness = 0.7, sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03, path_data = 'data')
skeletal = elem.generate(Revel)

cq.exporters.export(skeletal, 'skeletal.step')
