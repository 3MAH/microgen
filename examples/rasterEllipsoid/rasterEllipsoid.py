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
listPhases = []

elem = BasicGeometry(101,'ellipsoid', 0., 0.5, 0.5, 0, 0, 0, [0.15, 0.31, 0.4], '')
elli = elem.generate(Revel)

cq.exporters.export(elli, 'ellipsoid.step')
raster = rasterShapeList([elli],Revel,[5,5,5])

compound = cq.Compound.makeCompound(raster[0])
cq.exporters.export(compound, 'compound.step')

mesh('compound.step', raster[1], size_mesh, 1)

