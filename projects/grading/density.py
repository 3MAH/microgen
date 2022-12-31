import numpy as np
import pyvista as pv

from microgen import pv_tpms


def linear_graded_density(x, y, z):
    hmin = 0.3
    hmax = 1.5
    length = x[-1] - x[0]
    return (hmax - hmin) * x / length + 0.5 * (hmin + hmax)

def circular_graded_density(x: float, y: float, z: float) -> float:
    hmin = 0.3
    hmax = 1.5
    radius = 10.
    return (hmax - hmin) * (x**2 + y**2) / radius**2 + hmin


cell_size = np.array([5, 5, 1])
repeat_cell = np.array([5, 5, 1])
resolution = 20

geometry = pv_tpms.Tpms(
    surface_function=pv_tpms.gyroid,
    offset=linear_graded_density,
    repeat_cell=repeat_cell,
    cell_size=cell_size,
    resolution=resolution
)

plotter = pv.Plotter()
plotter.add_mesh(geometry.sheet, color="w")
plotter.camera_position = "xy"
plotter.parallel_projection = True
plotter.show_axes()
plotter.show()
