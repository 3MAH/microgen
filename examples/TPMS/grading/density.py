import numpy as np
import pyvista as pv

from microgen import Tpms
from microgen.shape.surface_functions import gyroid

def linear_graded_density(x, y, z):
    min_offset = 0.3
    max_offset = 1.5
    length = x[-1] - x[0]
    return (max_offset - min_offset) * x / length + 0.5 * (min_offset + max_offset)

def circular_graded_density(x: float, y: float, z: float) -> float:
    min_offset = 0.3
    max_offset = 2.0
    radius = 10.
    return (max_offset - min_offset) * (x**2 + y**2) / radius**2 + min_offset


geometry = Tpms(
    surface_function=gyroid,
    offset=circular_graded_density,
    repeat_cell=(5, 5, 1),
    cell_size=(5, 5, 1),
    resolution=30,
)
sheet = geometry.sheet

# plotter = pv.Plotter()
# plotter.add_mesh(geometry.sheet, color="w")
# plotter.view_xy()
# plotter.enable_parallel_projection()
# plotter.show_axes()
# plotter.show()
