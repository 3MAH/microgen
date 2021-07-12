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
    elem = phase(number[i],shape[i],xc[i],yc[i],zc[i],psi[i],theta[i],phi[i],param_geom,'')
    elem.generate(Revel)
    phases.append(elem)

#if periodicity==1:
#    for phase_elem in phases:
#        periodic(Revel,phase_elem)
#else:
#    for phase_elem in phases:
#        non_periodic(Revel,phase_elem)

for phase_elem in phases:
    p_phase = periodic(phase_elem,Revel)
    phases_perio.append(p_phase)

#solids = []
#for phase_list in phases_perio:
#    for phase in phase_list:
#        solids.append(phase.val())

sections = []
for phase_list in phases_perio:
    solids = []
    for phase in phase_list:
        solids.append(phase.val())
    section = fuse_parts(solids, False)
    sections.append(section)

#sections_cut = cut_parts(sections)
#print(sections_cut)
flat_list = [phase for phase_list in phases_perio for phase in phase_list]
print(len(flat_list))
#final_solid = fuse_parts(sections_cut, True)
final_solid = fuse_parts(sections, False)

#compound = cq.Compound.makeCompound(sections)
#compound = cq.Compound.makeCompound(sections_cut)

cq.exporters.export(final_solid, 'compound.step')
#cq.exporters.export(compound, 'compound.step')

#cq.Compound.exportBrep(final_solid, 'compound.brep')
#cq.Compound.exportBrep(compound, 'compound.brep')


#Mesh('compound.step', 0.1, 2)
MeshPeriodic('compound.step', 0.04, 2)

