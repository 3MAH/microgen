"""
Repeat a unit cell
==================

Import an octet-truss unit cell from a STEP file and repeat it on a
5 x 5 x 5 grid with :func:`microgen.repeat_shape`.

Requires the CAD backend: ``pip install microgen[cad]``.
"""

import pyvista as pv

from microgen import Rve, repeat_shape
from microgen.cad import import_step

rve = Rve(dim=1)

step_file = "octet_truss_cell.step"
unit_geom = import_step(step_file)

shape = repeat_shape(unit_geom, rve, grid=(5, 5, 5))

stl_file = "repeat_shape.stl"
shape.export_stl(stl_file)

# %%
pv.read(stl_file).plot(smooth_shading=True, split_sharp_edges=True, color="#ef7d2d")
