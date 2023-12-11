import numpy as np
import pyvista as pv

from microgen import Tpms
from microgen.shape.surface_functions import gyroid


def sigmoid(x, start, end):
    k = -0.4  # transition
    return start + (end - start) / (1 + np.exp(k * x))


def lerp(x, start, end):
    t = (x - np.min(x)) / (np.max(x) - np.min(x))
    return start + (end - start) * t


def gaussian(x, start, end):
    k = -0.02
    return start + (end - start) * np.exp(k * x**2)


def graded(x: float, y: float, z: float) -> float:
    min_cell_size = 2
    max_cell_size = 4
    # dim_x = dim_y = dim_z = lerp(x, min_cell_size, max_cell_size)
    dim_x = dim_y = dim_z = sigmoid(x, min_cell_size, max_cell_size)
    # dim_x = dim_y = dim_z = gaussian(x, max_cell_size, min_cell_size)
    return gyroid(
        x * max_cell_size / dim_x, y * max_cell_size / dim_y, z * max_cell_size / dim_z
    )


geometry = Tpms(
    surface_function=graded,
    offset=0.3,
    repeat_cell=(5, 5, 1),
    resolution=30,
)
sheet = geometry.sheet

# plotter = pv.Plotter()
# plotter.add_mesh(geometry.sheet, color="w")
# plotter.view_xy()
# plotter.enable_parallel_projection()
# plotter.show_axes()
# plotter.show()
