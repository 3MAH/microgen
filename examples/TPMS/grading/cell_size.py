import numpy as np
import pyvista as pv

from microgen import Tpms
from microgen.shape.surface_functions import gyroid


def sigmoid(x, start, end):
    k = -2  # transition
    return start + (end - start) / (1 + np.exp(k * x))


def lerp(x, start, end):
    t = (x - np.min(x)) / (np.max(x) - np.min(x))
    return start + (end - start) * t


def gaussian(x, start, end):
    k = -0.25
    return start + (end - start) * np.exp(k * x**2)


def graded(x: float, y: float, z: float) -> float:
    min_cell_size = 2
    max_cell_size = 4
    # qx = qy = qz = lerp(x, min_cell_size, max_cell_size)
    qx = qy = qz = sigmoid(x, min_cell_size, max_cell_size)
    # qx = qy = qz = gaussian(x, max_cell_size, min_cell_size)
    return gyroid(5 * x, qy * y, z)


geometry = Tpms(surface_function=graded, offset=0.3, cell_size=(5, 2, 1), resolution=100)
sheet = geometry.sheet

# plotter = pv.Plotter()
# plotter.add_mesh(geometry.sheet, color="w")
# plotter.view_xy()
# plotter.enable_parallel_projection()
# plotter.show_axes()
# plotter.show()
