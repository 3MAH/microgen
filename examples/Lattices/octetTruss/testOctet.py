import os
import numpy as np
import cadquery as cq
from microgen import removeEmptyLines, Rve, Cylinder, periodic, cutParts, meshPeriodic, Phase

# ----------LOADTXT------------------------------------------------------------------------------------------#

dir = os.path.dirname(os.path.realpath('__file__'))
# path
path_data = dir + '/'
Ngeomphase_file = 'test_octet.dat'

# fichier
NPhases_file = path_data + Ngeomphase_file
removeEmptyLines(NPhases_file)

dt = np.dtype([('number', int), ('shape', np.str_, 10),
               ('xc', np.float64), ('yc', np.float64), ('zc', np.float64),
               ('psi', np.float64), ('theta', np.float64), ('phi', np.float64),
               ('a1', np.float64), ('a2', np.float64)])
# précision du type des données
DATA = np.loadtxt(NPhases_file, dtype=dt,
                  usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                  skiprows=1, unpack=True, ndmin=1)

xc = DATA[2]
yc = DATA[3]
zc = DATA[4]
psi = DATA[5]
theta = DATA[6]
phi = DATA[7]
height = DATA[8]
radius = DATA[9]

# sections = read_sections(path_data,section_file)

rve = Rve(dim_x=1, dim_y=1, dim_z=1)
listPhases = []
listPeriodicPhases = []
n = len(xc)

for i in range(0, n):
    elem = Cylinder(center=(xc[i], yc[i], zc[i]),
                    orientation=(psi[i], theta[i], phi[i]),
                    height=height[i],
                    radius=radius[i])
    listPhases.append(Phase(shape=elem.generate()))

for phase_elem in listPhases:
    periodicPhase = periodic(phase=phase_elem, rve=rve)
    listPeriodicPhases.append(periodicPhase)

phases_cut = cutParts(cqShapeList=[phase.shape for phase in listPeriodicPhases], reverseOrder=False)
compound = cq.Compound.makeCompound([phase.shape for phase in phases_cut])

cq.exporters.export(compound, 'octettruss.step')
cq.exporters.export(compound, 'octettruss.stl')
meshPeriodic(mesh_file='octettruss.step', rve=rve, listPhases=phases_cut, size=0.03, order=1, output_file='octettruss.vtk')
