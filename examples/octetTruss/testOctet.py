import os
import sys
import numpy as np
import cadquery as cq
import gmsh
from microgen import *

#----------LOADTXT------------------------------------------------------------------------------------------#

dir = os.path.dirname(os.path.realpath('__file__'))
#path
path_data = dir + '/'
Ngeomphase_file = 'test_octet.dat'

#fichier
NPhases_file = path_data + Ngeomphase_file
removeEmptyLines(NPhases_file)

dt = np.dtype([('number',int),('shape', np.str_, 10), ('xc', np.float64),('yc', np.float64),('zc', np.float64),
               ('psi', np.float64),('theta', np.float64),('phi', np.float64),('a1', np.float64),('a2', np.float64)])
#précision du type des données
number,shape,xc,yc,zc,psi,theta,phi,a1,a2,= np.loadtxt(NPhases_file,dtype=dt, usecols=(0,1,2,3,4,5,6,7,8,9),
                                             skiprows=1,unpack=True, ndmin=1)
   
#sections = read_sections(path_data,section_file)

#Size of the mesh
size_mesh = 0.03
a = 1.0
b = 1.0
c = 1.0

Revel = Rve(a,b,c,size_mesh)
listPhases = []
listPeriodicPhases = []
n=len(xc)

for i in range (0,n):
    param_geom=[a1[i],a2[i]]
    elem = BasicGeometry(number[i],shape[i],xc[i],yc[i],zc[i],psi[i],theta[i],phi[i],param_geom,'')
    listPhases.append(elem.generate())

for phase_elem in listPhases:
    print(phase_elem)
    periodicPhase = periodic(phase_elem,Revel)
    listPeriodicPhases.append(periodicPhase)

phases_cut = cutParts([s[0] for s in listPeriodicPhases], False)
compound = cq.Compound.makeCompound(phases_cut[0])

cq.exporters.export(compound, 'compound.step')
#mesh('compound.step', phases_cut[1], 0.03, 1)
meshPeriodic('compound.step', Revel, phases_cut[1], 0.03, 1)
