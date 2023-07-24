import numpy as np
import pyvista as pv

from microgen import CylindricalTpms
from microgen.shape.surface_functions import gyroid


def swapped_gyroid(x, y, z):
    return gyroid(x=z, y=y, z=x)

def grading(x, y, z):
    a = (max_offset - min_offset) / (np.max(x) - np.min(x))
    b = max_offset - a * np.max(x)

    return a * x + b


min_offset = 0.
max_offset = 3.

geometry = CylindricalTpms(
    radius=10,
    surface_function=swapped_gyroid,
    offset=grading,
    cell_size=(1, 1, 1),
    repeat_cell=(10, 0, 1),
    resolution=20,
)
sheet = geometry.sheet

# plotter = pv.Plotter()
# plotter.add_mesh(sheet, color="w")
# plotter.view_xy()
# plotter.enable_parallel_projection()
# plotter.show_axes()
# plotter.show()
