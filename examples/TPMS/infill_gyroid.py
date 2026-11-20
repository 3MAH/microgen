"""
Gyroid infilled by a gyroid
===========================

A gyroid sheet used as the envelope of a finer gyroid infill: a hierarchical
TPMS.
"""

import pyvista as pv

from microgen import Infill, Tpms
from microgen.shape.surface_functions import gyroid

tpms = Tpms(
    surface_function=gyroid,
    offset=1.0,
    resolution=30,
)

infill = Infill(
    obj=tpms.sheet,
    surface_function=gyroid,
    cell_size=0.1,
    offset=0.5,
    resolution=15,
)

# %%
pl = pv.Plotter()
pl.add_mesh(infill.sheet, color="#ef7d2d")
pl.add_axes()
pl.show()
