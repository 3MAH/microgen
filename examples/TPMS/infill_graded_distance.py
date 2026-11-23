"""
Infill graded by the distance to the surface
============================================

Fill the Stanford bunny with a gyroid whose offset is given by
:class:`microgen.NormedDistance`: thin near the skin, thick in the core. The
infill is clipped to show the inside.
"""

import pyvista as pv
from pyvista import examples

from microgen import Infill, NormedDistance
from microgen.shape.surface_functions import gyroid

mesh = examples.download_bunny()

offset = NormedDistance(
    obj=mesh,
    boundary_offset=0,
    furthest_offset=3.0,
    boundary_weight=1.0,
)

infill = Infill(
    obj=mesh,
    surface_function=gyroid,
    cell_size=0.015,
    offset=offset,
    resolution=5,
)

# %%
pl = pv.Plotter()
pl.add_mesh(infill.sheet.clip("z", origin=mesh.center_of_mass()), color="#ef7d2d")
pl.add_axes()
pl.view_xy()
pl.show()
