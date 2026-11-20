"""
Spherical gyroid shell
======================

A gyroid surface mapped onto a spherical coordinate frame with
:class:`microgen.SphericalTpms`.
"""

import pyvista as pv

from microgen import SphericalTpms
from microgen.shape.surface_functions import gyroid


def swapped_gyroid(x, y, z):
    return gyroid(x=z, y=y, z=x)


geometry = SphericalTpms(
    surface_function=swapped_gyroid,
    offset=0.5,
    repeat_cell=(1, 10, 20),
    radius=3,
    resolution=50,
)
surface = geometry.surface

# %%
plotter = pv.Plotter()
plotter.add_mesh(surface, smooth_shading=True, color="#ef7d2d")
plotter.view_xy()
plotter.enable_parallel_projection()
plotter.show_axes()
plotter.show()
