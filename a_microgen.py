# -*- coding: utf8 -*-
import os
import sys
import numpy as np
import cadquery as cq
import gmsh
from microgen import *


#----------LOADTXT------------------------------------------------------------------------------------------#

dir = os.path.dirname(os.path.realpath('__file__'))
#chemin
path_data = dir + '/'
#Ngeomphase_file = 'test_scaffold.dat'
#section_file = 'Nsections_scaffold.dat'
#Ngeomphase_file = 'test_diamond.dat'
#section_file = 'Nsections_diamond.dat'
#Ngeomphase_file = 'test_scaffbrid.dat'
#section_file = 'Nsections_scaffbrid.dat'
#Ngeomphase_file = 'test_electrospin.dat'
#section_file = 'Nsections_electrospin.dat'
#Ngeomphase_file = 'test_electrospin_mat.dat'
#section_file = 'Nsections_electrospin_mat.dat'
#Ngeomphase_file = 'test_spheres.dat'
#section_file = 'Nsections_spheres.dat'
Ngeomphase_file = 'test_octet.dat'
#Ngeomphase_file = 'test_octet_slim.dat'
section_file = 'Nsections_octet.dat'
#Ngeomphase_file = 'test_octet_spheres.dat'
#section_file = 'Nsections_octet_spheres.dat'
#Ngeomphase_file = 'test_fibers.dat'
#section_file = 'Nsections_fibers.dat'
#Ngeomphase_file = 'test_sphere.dat'
#section_file = 'Nsections.dat'

#fichier
NPhases_file = path_data + Ngeomphase_file
remove_empty_lines(NPhases_file)

dt = np.dtype([('number',int),('shape', np.str_, 10), ('xc', np.float64),('yc', np.float64),('zc', np.float64),
               ('psi', np.float64),('theta', np.float64),('phi', np.float64),('a1', np.float64),('a2', np.float64)])
#précision du type des données
number,shape,xc,yc,zc,psi,theta,phi,a1,a2,= np.loadtxt(NPhases_file,dtype=dt, usecols=(0,1,2,3,4,5,6,7,8,9),
                                             skiprows=1,unpack=True, ndmin=1)
   
sections = read_sections(path_data,section_file)

#Size of the mesh
size_mesh = 0.03
a = 1.0
b = 1.0
c = 1.0

periodicity = 0

Revel = rve(a,b,c,size_mesh)
phases = []
phases_perio = []
n=len(xc)
for i in range (0,n):
    param_geom=[a1[i],a2[i]]
    elem = BasicGeometry(number[i],shape[i],xc[i],yc[i],zc[i],psi[i],theta[i],phi[i],param_geom,'')
    phases.append(elem.generate(Revel))

for phase_elem in phases:
    print(phase_elem)
    p_phase = periodic(phase_elem,Revel)
    print(p_phase)
    phases_perio.append(p_phase)

print(phases_perio)

phases_cut = cut_parts([s[0] for s in phases_perio], False)
compound = cq.Compound.makeCompound(phases_cut[0])

#phases_cut = cut_parts(phases)
#compound = cq.Compound.makeCompound(phases_cut[0])

#cq.exporters.export(final_solid[0], 'compound.step')

#cq.Compound.exportBrep(final_solid, 'compound.brep')
#cq.Compound.exportBrep(compound, 'compound.brep')

cq.exporters.export(compound, 'compound.step')
#Mesh('compound.step', phases_cut[1], 0.03, 1)
MeshPeriodic('compound.step', phases_cut[1], 0.03, 1)

#cq.exporters.export(compound, 'compound.step')
#Mesh('compound.step', [s[1] for s in phases_perio], 0.05, 1)


#MeshPeriodic('compound.step', 0.04, 2)

