import numpy as np
import pyvista as pv

from microgen import Tpms
from microgen.shape.surface_functions import gyroid, schwarzP, fischerKochS

repeat = 2
k = 5
point_1 = (-0.5, 0, 0)
point_2 = (0.5, 0, 0)
point_3 = (0, 0.5, 0)


def exp_func(x, y, z):
    norm = x**2 + y**2 + z**2
    return 1 + np.exp(k * norm)

def weight(x, y, z, points, index):
    denom = 0.0
    for p in points:
        denom += exp_func(x - p[0], y - p[1], z - p[2])
    point = points[index]
    return exp_func(x - point[0], y - point[1], z - point[2]) / denom

def multi_morph(phi, x, y, z, points):
    result = 0.0
    for index, surface_function in enumerate(phi):
        weight_func = weight(x, y, z, points, index)
        result += weight_func * surface_function(repeat * x, repeat * y, repeat * z)
    return result

def trigraded(x, y, z):
    return multi_morph(
        phi=[schwarzP, gyroid, fischerKochS],
        x=x, y=y, z=z,
        points=[point_1, point_2, point_3]
    )


geometry = Tpms(
    surface_function=trigraded,
    offset=0.3,
    cell_size=1.,
    repeat_cell=repeat,
    resolution=100
)

geometry.sheet.save("trigraded.vtk")
# plotter = pv.Plotter()
# plotter.add_mesh(geometry.sheet, color="w")
# plotter.view_xy()
# plotter.enable_parallel_projection()
# plotter.show_axes()
# plotter.show()
