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

elem = BasicGeometry(0,'tpms',0.5,0.5,0.5,0.,0.,0.,['gyroid','sheet','na'],'data')
#elem = BasicGeometry(0,'tpms',0.5,0.5,0.5,0.,0.,0.,['gyroid','skeletal','na'],'data')

elem.geometry.createSurfaces('data', Revel, thickness = 0.3, sizeMesh=0.05, minFacetAngle=10., maxRadius=0.05)
elem.createGeometry(Revel)

temp = elem.generate(Revel)
phases.append(temp)

#temp2 = temp.translate((0, 0, 1))
#phases.append(temp2)

occ_solids_list = [s.Solids() for s in phases]
flat_list = [item.copy() for sublist in occ_solids_list for item in sublist]
print(flat_list)

compound = cq.Compound.makeCompound(phases)
cq.exporters.export(compound, 'compound.step')

Mesh('compound.step', [flat_list,],  0.03, 1)
#MeshPeriodic('compound.step', Revel, [flat_list,],  0.04, 1)
