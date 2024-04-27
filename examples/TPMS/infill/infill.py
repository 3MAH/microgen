"""TPMS infill using the Infill class."""

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

pl = pv.Plotter()
pl.add_mesh(infill.sheet, color="w")
pl.add_axes()
# pl.show()
