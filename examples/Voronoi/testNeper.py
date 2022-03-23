import os
import sys
import numpy as np
import cadquery as cq
import gmsh
from microgen import *

#We import the Polyhedra from Nepoer tesselation file
listPolyhedra, seed, vertices, edges, faces, polys = parseNeper('test1')

#We fix the size of the mesh and dimensions of the RVE
size_mesh = 0.05
a = 1.0
b = 1.0
c = 1.0

periodicity = 0

Revel = Rve(a,b,c,size_mesh)
phases = []
listSolids = []

for Polyhedra in listPolyhedra:
    elem = BasicGeometry(0,'polyhedra',Polyhedra["original"][0],Polyhedra["original"][1], Polyhedra["original"][2],0.,0.,0.,Polyhedra,'')
    temp = elem.generate()
    phases.append(temp)

compound = cq.Compound.makeCompound(phases)
cq.exporters.export(compound, 'compound.step')

occ_solids_list = [s.Solids() for s in phases]
mesh('compound.step', occ_solids_list, size_mesh, 1)
