import microgen
import cadquery as cq
import numpy as np
import os

def test_octettruss():

    # fichier
    NPhases_file = "examples/octetTruss/test_octet.dat"
    microgen.removeEmptyLines(NPhases_file)

    dt = np.dtype([('number', int), ('shape', np.str_, 10),
                ('xc', np.float64), ('yc', np.float64), ('zc', np.float64),
                ('psi', np.float64), ('theta', np.float64), ('phi', np.float64),
                ('height', np.float64), ('radius', np.float64)])
    # précision du type des données
    number, shape, xc, yc, zc, psi, theta, phi, height, radius, = np.loadtxt(NPhases_file, dtype=dt,
                                                                             usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                                                                             skiprows=1, unpack=True, ndmin=1)

    # Size of the mesh
    size_mesh = 0.03
    a = 1.0
    b = 1.0
    c = 1.0

    revel = microgen.Rve(a, b, c, size_mesh)
    listPhases = []
    listPeriodicPhases = []
    n = len(xc)

    for i in range(0, n):
        elem = microgen.BasicGeometry(number=number[i], shape=shape[i],
                            xc=xc[i], yc=yc[i], zc=zc[i],
                            psi=psi[i], theta=theta[i], phi=phi[i],
                            param_geom={"height": height[i],
                                        "radius": radius[i]},
                            path_data='')
        listPhases.append(elem.generate())

    for phase_elem in listPhases:
        print(phase_elem)
        periodicPhase = microgen.periodic(cqshape=phase_elem, rve=revel)
        listPeriodicPhases.append(periodicPhase)

    phases_cut = microgen.cutParts(cqShapeList=[s[0] for s in listPeriodicPhases], reverseOrder=False)
    compound = cq.Compound.makeCompound(phases_cut[0])

    os.makedirs('tests/data', exist_ok=True) # if data folder doesn't exist yet

    cq.exporters.export(compound, 'tests/data/compound.step')
    microgen.meshPeriodic(mesh_file='tests/data/compound.step', rve=revel, listPhases=phases_cut[1], size=0.03, order=1)