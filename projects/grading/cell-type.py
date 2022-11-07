from microgen import Tpms, tpms

import numpy as np
import pyvista as pv
import cadquery as cq


# def G(x, y, z):
#     t = 2 # radius
#     return (x - 2)**2 + y**2 + z**2 - t**2
#     # return x * (5 - x)
#     # return x**2 - t**2

# def gamma(x, y, z):
#     k = -2 # transition
#     return 1. / (1 + np.exp(k*G(x, y, z)))

# def graded(x: float, y: float, z: float, q: tuple[int, int, int], l: tuple[float, float, float]) -> float:
#     k = 2
#     ll = (k*l[0], k*l[1], k*l[2])
#     # return gamma(x, y, z) * gyroid(x, y, z, q, l) + (1 - gamma(x, y, z)) * honeycomb(x, y, z, q, l)
#     return gamma(x, y, z) * tpms.gyroid(x, y, z, q, l) + (1 - gamma(x, y, z)) * tpms.gyroid(x, y, z, q, ll)
#     # L = l[0]
#     # y1 = gyroid(x, y, z, q, l)
#     # y2 = gyroid(x, y, z, q, ll)
#     # return (y2 - y1) * x / L + 0.5 * (y1 + y2)


# geometry = Tpms(
#     surface_function=graded,
#     type_part="sheet",
#     thickness=0.2,
#     repeat_cell=(8, 16, 16),
#     cell_size=(2, 4, 4)
# )
# # shape = geometry.cylinder(isovalue=0., resolution=20, smoothing=0)
# # shape = geometry.generate(resolution=100, smoothing=0)
# shape = geometry.generateSurfaceVtk(isovalue=0, resolution=400, smoothing=0)
# p = pv.Plotter()
# p.add_mesh(shape, color='w')
# p.camera_position = 'xy'
# p.parallel_projection = True 
# p.show()

# # geometry = Tpms(
# #     surface_function=tpms.gyroid,
# #     type_part=tpms.SHEET,
# #     thickness=0.2
# # )
# # shape = geometry.generate(resolution=20, smoothing=0)
# # cq.exporters.export(shape, "test.stl")


l = 1
k = 5
p1 = (-0.5, 0, 0)
p2 = (0.5, 0, 0)
p3 = (0, 0.5, 0)
phi1 = tpms.schwarzP
phi2 = tpms.gyroid
def phi3(x, y, z, q, l):
    kx = 2* np.pi * q[0] / l[0]
    ky = 2* np.pi * q[1] / l[1]
    kz = 2* np.pi * q[2] / l[2]
    return np.cos(kx * x) * np.cos(ky * y) * np.cos(kz * z) - np.sin(kx * x) * np.sin(ky * y) * np.sin(kz * z)

def exp_func(x, y, z):
    norm = x**2 + y**2 + z**2
    return 1 + np.exp(k * norm)

def weight(x, y, z, p, i):
    denom = 0.
    for j in range(len(p)):
        denom += exp_func(x - p[j][0], y - p[j][1], z - p[j][2])
    return exp_func(x - p[i][0], y - p[i][1], z - p[i][2]) / denom

def multi_morph(phi, x, y, z, p, q, l):
    result = 0.
    for i in range(len(phi)):
        weight_func = weight(x, y, z, p, i)
        result += weight_func * phi[i](x, y, z, q, l)
    return result

def trigraded(x, y, z, q, l):
    return multi_morph(phi=[phi1, phi2, phi3], x=x, y=y, z=z, p=[p1, p2, p3], q=q, l=l)

geometry = Tpms(
    surface_function=tpms.gyroid,
    type_part="sheet",
    thickness=0.1,
    repeat_cell=5,
    cell_size=5
)
# shape = geometry.cylinder(isovalue=0., resolution=20, smoothing=0)
# shape = geometry.generate(resolution=100, smoothing=0)
shape = geometry.generate(resolution=100)
cq.exporters.export(shape, "test.stl")
# p = pv.Plotter()
# p.add_mesh(shape, color='w')
# p.camera_position = 'xy'
# p.parallel_projection = True 
# p.show_axes()
# p.show()