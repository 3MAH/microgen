"""
Periodic octet-truss mesh
=========================

Build an octet-truss from cylinders read in a data file, make it periodic
with :func:`microgen.periodic_split_and_translate`, and mesh it with
:func:`microgen.mesh_periodic` so that opposite faces have matching nodes.

Requires the CAD backend: ``pip install microgen[cad]``.
"""

import numpy as np
import pyvista as pv

from microgen import (
    Cylinder,
    Phase,
    Rve,
    cut_phases,
    mesh_periodic,
    periodic_split_and_translate,
)
from microgen.cad import make_compound

# %%
# Read the strut positions, orientations and dimensions.
data_file = "octet_truss.dat"

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
        ("a1", np.float64),
        ("a2", np.float64),
    ]
)
data = np.loadtxt(
    data_file,
    dtype=dt,
    usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
    skiprows=1,
    unpack=True,
    ndmin=1,
)

xc = data[2]
yc = data[3]
zc = data[4]
psi = data[5]
theta = data[6]
phi = data[7]
height = data[8]
radius = data[9]

# %%
# Build the struts, split the parts that cross the RVE faces and translate
# them to the opposite side.
rve = Rve(dim=1)
list_phases = []
for i in range(len(xc)):
    elem = Cylinder(
        center=(xc[i] - 0.5, yc[i] - 0.5, zc[i] - 0.5),
        orientation=(psi[i], theta[i], phi[i]),
        height=height[i],
        radius=radius[i],
    )
    list_phases.append(Phase.from_cad(elem.generate_cad()))

periodic_phases = [
    periodic_split_and_translate(phase=phase, rve=rve) for phase in list_phases
]

phases_cut = cut_phases(phases=periodic_phases, reverse_order=False)
compound = make_compound([phase.cad for phase in phases_cut])

step_file = "octet_truss.step"
compound.export_step(step_file)
compound.export_stl("octet_truss.stl")

# %%
# Periodic tetrahedral mesh.
vtk_file = "octet_truss.vtk"
mesh_periodic(
    mesh_file=step_file,
    rve=rve,
    list_phases=phases_cut,
    order=1,
    size=0.03,
    output_file=vtk_file,
)

# %%
pv.read(vtk_file).plot(show_edges=True, color="#ef7d2d")
