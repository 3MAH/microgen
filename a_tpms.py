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

Revel = rve(a,b,c,size_mesh)
phases = []

elem = BasicGeometry(0,'tpms',0.5,0.5,0.5,0.,0.,0.,['gyroid','sheet','na'],'')
temp = elem.generate(Revel)
phases.append(temp)

compound = cq.Compound.makeCompound(phases)
cq.exporters.export(compound, 'compound.step')
#Mesh('compound.step', phases_cut[1], 0.03, 1)
