from __future__ import annotations
from functools import partial

import numpy as np
import pyvista as pv

from microgen import CylindricalTpms
from microgen.shape.surface_functions import gyroid


def swapped_gyroid(x, y, z):
    return gyroid(x=z, y=y, z=x)


def grading(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    min_offset: float | np.ndarray = 0.0,
    max_offset: float | np.ndarray = 3.0,
) -> np.ndarray:
    rad = np.sqrt(x**2 + y**2)
    max_rad = np.max(rad)
    min_rad = np.min(rad)

    return min_offset + (max_offset - min_offset) * (rad - min_rad) / (
        max_rad - min_rad
    )


full_density_offset = CylindricalTpms.offset_from_density(
    surface_function=swapped_gyroid,
    part_type="sheet",
    density=1.0,
)

geometry = CylindricalTpms(
    radius=5,
    surface_function=swapped_gyroid,
    offset=partial(
        grading,
        min_offset=0.0,
        max_offset=full_density_offset,
    ),
    cell_size=(1, 1, 1),
    repeat_cell=(5, 0, 1),
    resolution=20,
)
sheet = geometry.sheet

# plotter = pv.Plotter()
# plotter.add_mesh(sheet, color="w")
# plotter.view_xy()
# plotter.enable_parallel_projection()
# plotter.show_axes()
# plotter.show()
