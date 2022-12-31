from typing import List

import numpy as np
import pyvista as pv


def swapped_gyroid(x, y, z):
    return np.cos(x) * np.sin(y) + np.cos(y) * np.sin(z) + np.cos(z) * np.sin(x)


repeat_cell = np.array([1, 1, 1])
cell_size = np.array([1, 1, 1])
center_offset = 0.5
resolution = 20

offset = 1

linspaces: List[np.ndarray] = [
    np.linspace(
        -center_offset * cell_size_axis * repeat_cell_axis,
        center_offset * cell_size_axis * repeat_cell_axis,
        resolution * repeat_cell_axis,
    )
    for repeat_cell_axis, cell_size_axis in zip(repeat_cell, cell_size)
]

x, y, z = np.meshgrid(*linspaces)

grid = pv.StructuredGrid(x, y, z)
kx, ky, kz = 2 * np.pi / cell_size

surface_function = swapped_gyroid(kx * x, ky * y, kz * z)

grid["lower_surface"] = (surface_function - offset / 2.0).ravel(order="F")
grid["upper_surface"] = (surface_function + offset / 2.0).ravel(order="F")
sheet = (
    grid.clip_scalar(scalars="upper_surface", invert=False)
    .clip_scalar(scalars="lower_surface")
    .extract_surface()
)
upper_skeletal = grid.clip_scalar(scalars="upper_surface")
lower_skeletal = grid.clip_scalar(scalars="lower_surface", invert=False)

print(f"relative density = {sheet.volume / grid.volume:.2%}")

pl = pv.Plotter()
pl.add_mesh(sheet, color="w")
# pl.add_mesh(upper_skeletal, color="b")
# pl.add_mesh(lower_skeletal, color="r")
pl.show()
