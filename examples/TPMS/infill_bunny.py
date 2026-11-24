"""
Gyroid infill of a mesh
=======================

Fill the Stanford bunny with a gyroid with :class:`microgen.Infill`.
"""

import pyvista as pv
from pyvista import examples

from microgen import Infill
from microgen.shape.surface_functions import gyroid

bunny = examples.download_bunny()

infill = Infill(
    obj=bunny,
    surface_function=gyroid,
    cell_size=0.015,
    offset=0.5,
    resolution=10,
)

# %%
pl = pv.Plotter()
pl.add_mesh(infill.sheet, color="#ef7d2d")
pl.add_axes()
pl.view_xy()
pl.show()
