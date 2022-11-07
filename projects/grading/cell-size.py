from microgen import Tpms, tpms

import numpy as np
import pyvista as pv


def sigmoid(x, q1, q2):
    k = -2 # transition
    return q1 + (q2 - q1) / (1 + np.exp(k*x))

def linear(x, q1, q2, lx):
    a = (q2 - q1) / lx
    b = (q1 + q2) * 0.5
    return a * x + b

def gaussian(x, q1, q2):
    k = -1.5
    return q1 + (q2 - q1) * np.exp(k * x**2)


def graded(x: float, y: float, z: float, q: tuple[int, int, int], l: tuple[float, float, float]) -> float:
    q1 = 20
    q2 = 10
    # qx = qy = qz = linear(x, q1, q2, l[0])
    # qx = qy = qz = sigmoid(x, q1, q2)
    qx = qy = qz = gaussian(x, q1, q2)
    qq = (qx, qy, qz)
    # return tpms.schwarzP(x, y, z, qq, l)
    return tpms.gyroid(x, y, z, qq, l)


geometry = Tpms(
    surface_function=graded,
    type_part="sheet",
    thickness=0.2,
    repeat_cell=(8, 8, 8),
    cell_size=(4, 4, 4)
)

shape = geometry.generateSurfaceVtk(isovalue=0, resolution=300, smoothing=0)
p = pv.Plotter()
p.add_mesh(shape, color='w')
p.camera_position = 'xy'
p.parallel_projection = True 
p.show_axes()
p.show()