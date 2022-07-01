import microgen
import cadquery as cq
import numpy as np
import os

import pytest


@pytest.mark.filterwarnings("ignore:Object intersecting")
def test_octettruss():

    # fichier
    NPhases_file = "examples/Lattices/octetTruss/test_octet.dat"

    dt = np.dtype(
        [
            ("number", int),
            ("shape", np.str_, 10),
            ("xc", np.float64),
            ("yc", np.float64),
            ("zc", np.float64),
            ("psi", np.float64),
            ("theta", np.float64),
            ("phi", np.float64),
            ("height", np.float64),
            ("radius", np.float64),
        ]
    )
    # précision du type des données
    number, shape, xc, yc, zc, psi, theta, phi, height, radius, = np.loadtxt(
        NPhases_file,
        dtype=dt,
        usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
        skiprows=1,
        unpack=True,
        ndmin=1,
    )

    revel = microgen.Rve(dim_x=1, dim_y=1, dim_z=1, center=(0.5, 0.5, 0.5))
    listPhases = []  # type: list[microgen.Phase]
    listPeriodicPhases = []  # type: list[microgen.Phase]
    n = len(xc)

    for i in range(0, n):
        elem = microgen.shape.cylinder.Cylinder(
            center=(xc[i], yc[i], zc[i]),
            orientation=(psi[i], theta[i], phi[i]),
            height=height[i],
            radius=radius[i],
        )
        phase = microgen.Phase(shape=elem.generate())
        listPhases.append(phase)

    for phase_elem in listPhases:
        periodicPhase = microgen.periodic(phase=phase_elem, rve=revel)
        listPeriodicPhases.append(periodicPhase)

    print([phase.centerOfMass for phase in listPeriodicPhases])

    phases_cut = microgen.cutPhases(phaseList=listPeriodicPhases, reverseOrder=True)
    phases_cut = microgen.cutPhases(phaseList=listPeriodicPhases, reverseOrder=False)
    compound = cq.Compound.makeCompound([phase.shape for phase in phases_cut])

    os.makedirs("tests/data", exist_ok=True)  # if data folder doesn't exist yet

    cq.exporters.export(compound, "tests/data/compound.step")
    microgen.meshPeriodic(
        mesh_file="tests/data/compound.step",
        rve=revel,
        listPhases=phases_cut,
        size=0.03,
        order=1,
        output_file="tests/data/MeshPeriodic.vtk",
    )


if __name__ == "__main__":
    test_octettruss()
