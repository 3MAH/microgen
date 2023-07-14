import numpy as np
import pyvista as pv

from microgen import CylindricalTpms
from microgen.shape.surface_functions import gyroid


def swapped_gyroid(x, y, z):
    return gyroid(x=y, y=z, z=x)

def grading(x, y, z):
    a = (max_offset - min_offset) / (np.max(x) - np.min(x))
    b = max_offset - a * np.max(x)

    return a * x + b

geometry = CylindricalTpms(
    radius=10,
    surface_function=swapped_gyroid,
    offset=grading,
    cell_size=(1, 1, 1),
    repeat_cell=(10, 0, 1),
)


min_offset = 0.
max_offset = 3.

plotter = pv.Plotter()
plotter.add_mesh(geometry.sheet, color="w")
plotter.view_xy()
plotter.enable_parallel_projection()
plotter.show_axes()
plotter.show()
