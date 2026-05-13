"""Gyroid infilled by gyroid."""

import pyvista as pv

from microgen import Infill, Tpms
from microgen.shape.surface_functions import gyroid

tpms = Tpms(surface_function=gyroid).with_offset(1.0).with_resolution(30)

infill = (
    Infill(obj=tpms.sheet, surface_function=gyroid)
    .with_cell_size(0.1)
    .with_offset(0.5)
    .with_resolution(15)
)

pl = pv.Plotter()
pl.add_mesh(infill.sheet, color="w")
pl.add_axes()
# pl.show()
