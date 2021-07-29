import os
import sys
import numpy as np
import cadquery as cq
import gmsh
from microgen import *
import pyvoro
from random import *


AL = []
for i in range(500) :
   L=[]
   for j in range(3) :
      a=random()
     # a=a*10
      L.append(a)
   AL.append(L)

A = pyvoro.compute_voronoi(
 # [[1.0, 2.0, 3.0], [4.0, 5.5, 6.0], [7.0, 4.5, 6.0],[9.0, 7.5, 9.0],[1.0, 1.5, 1.0],[1.0, 1.5, 9.0],[1.0, 9.5, 9.0],[9.0, 1.5, 9.0],[5.0, 1.5, 5.0],[1.0, 5.5, 1.0],[2.0, 1.0, 8.0], [6.0, 1.5, 9.0], [4.0, 5.5, 2.0],[7.0, 1.5, 9.0],[0.5, 1.5, 6.0],[0.5, 2.5, 8.0],[1.0, 5.5, 9.0],[9.0, 9.5, 7.0],[5.0, 9.5, 2.0],[1.0, 6.5, 3.0]], # point positions
   AL,
  [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], # limits
  1.0, # block size
  #radii=[1.3, 1.4] # particle radii -- optional, and keyword-compatible arg.
)

#Size of the mesh
size_mesh = 0.03
a = 1.0
b = 1.0
c = 1.0

periodicity = 0

Revel = rve(a,b,c,size_mesh)
phases = []

for voro in A:
    elem = BasicGeometry(0,'voronoi',voro["original"][0],voro["original"][1], voro["original"][2],0.,0.,0.,voro,'')
    temp = elem.generate(Revel)
#    c0 = temp.Center()
#    transform_mat = cq.Matrix([[0.9, 0, 0, c0.x],
#    [0, 0.9, 0, c0.y],
#    [0, 0, 0.9, c0.z]])
#    temp2 = temp.transformGeometry(transform_mat)
#    c1 = temp2.Center()
    phases.append(temp)

compound = cq.Compound.makeCompound(phases)
cq.exporters.export(compound, 'compound.step')
#Mesh('compound.step', phases_cut[1], 0.03, 1)
