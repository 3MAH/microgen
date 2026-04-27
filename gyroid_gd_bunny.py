from typing import List

import numpy as np
import pyvista as pv
from pyvista import examples


def gyroid(x, y, z):
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


repeat_cell = np.array([8, 8, 8])
cell_size = np.array([1, 1, 1])
center_offset = 0.5
resolution = 15

offset = np.pi / 2.0

linspaces: List[np.ndarray] = []
for repeat_cell_axis, cell_size_axis in zip(repeat_cell, cell_size):
    linspaces.append(
        np.linspace(
            -center_offset * cell_size_axis * repeat_cell_axis,
            center_offset * cell_size_axis * repeat_cell_axis,
            resolution * repeat_cell_axis,
        )
    )

x, y, z = np.meshgrid(*linspaces)

grid = pv.StructuredGrid(x, y, z)
kx, ky, kz = 2 * np.pi / cell_size

surface_function = gyroid(kx * x, ky * y, kz * z)

bunny = examples.download_bunny()
transform_matrix = np.array([[40, 0, 0, 0], [0, 40, 0, 0], [0, 0, 40, 0], [0, 0, 0, 1]])
bunny.transform(transform_matrix, inplace=True)
center = bunny.center_of_mass()
bunny.translate(-center, inplace=True)
print("newcenter = ", bunny.center_of_mass())
print(bunny.bounds)
print(grid.bounds)

grid.compute_implicit_distance(bunny, inplace=True)

# normalize :
dist = -1.0 * grid["implicit_distance"]
dist[dist > 0] = 0
dist_norm = (dist - min(dist)) / (max(dist) - min(dist))
x_t = 0.5
l = 0.2
reg_func = 0.6 * (1.0 + np.tanh((dist_norm - x_t) / l)) - 0.2

print(min(dist))
print(max(dist))

print(min(dist_norm))
print(max(dist_norm))

print(min(reg_func))
print(max(reg_func))

grid["lower_surface"] = surface_function.ravel(order="F") - offset * reg_func
grid["upper_surface"] = surface_function.ravel(order="F") + offset * reg_func
sheet = grid.clip_scalar(scalars="upper_surface", invert=False).clip_scalar(
    scalars="lower_surface"
)
upper_skeletal = grid.clip_scalar(scalars="upper_surface")
lower_skeletal = grid.clip_scalar(scalars="lower_surface", invert=False)

clipped = sheet.clip_surface(bunny, invert=False)
clipped.compute_implicit_distance(bunny, inplace=True)

print(f"relative density = {clipped.volume / bunny.volume:.2%}")

clipped2 = clipped.clip("y", origin=(0, -0.5, 0.0), invert=False)
clipped2["dist"] = -1.0 * clipped2["implicit_distance"]

pl = pv.Plotter()
# pl.add_mesh(grid, color="b", opacity = 0.1)
pl.add_mesh(bunny, color="w", opacity=0.1)
pl.add_mesh(clipped2, scalars="dist", cmap="inferno")
pl.show()
