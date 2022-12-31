import numpy as np
import pyvista as pv

from microgen import Tpms, pv_tpms


def sigmoid(x, q1, q2):
    k = -2  # transition
    return q1 + (q2 - q1) / (1 + np.exp(k * x))


def lerp(x, q1, q2):
    t = (x - min(x)) / (max(x) - min(x))
    return q1 + (q2 - q1) * t


def gaussian(x, q1, q2):
    k = -1.5
    return q1 + (q2 - q1) * np.exp(k * x**2)


def graded(
    x: float, y: float, z: float) -> float:
    q1 = 6
    q2 = 3
    # qx = qy = qz = lerp(x, q1, q2)
    qx = qy = qz = sigmoid(x, q1, q2)
    # qx = qy = qz = gaussian(x, q1, q2)
    # return tpms.schwarzP(qx * x, qy * y, qz * z)
    return pv_tpms.gyroid(qx * x, qy * y, z)


# geometry = Tpms(surface_function=graded)
geometry = Tpms(surface_function=graded, cell_size=(6, 6, 1), resolution=100)

p = pv.Plotter()
p.add_mesh(geometry.surface, color="w")
p.camera_position = "xy"
p.parallel_projection = True
p.show_axes()
p.show()
