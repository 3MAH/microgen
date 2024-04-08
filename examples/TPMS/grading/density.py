"""Example of density grading of a TPMS geometry.

**Note:** The density does not vary linearly with the offset.
"""

import numpy as np
import pyvista as pv

from microgen import Tpms
from microgen.shape.surface_functions import gyroid


def linear_graded_offset(x: np.ndarray, _: np.ndarray, __: np.ndarray) -> np.ndarray:
    """Linearly graded offset."""
    min_offset = 0.0
    max_offset = 3.0
    length = np.max(x) - np.min(x)
    return (max_offset - min_offset) * x / length + 0.5 * (min_offset + max_offset)


def circular_graded_offset(x: np.ndarray, y: np.ndarray, _: np.ndarray) -> np.ndarray:
    """Circularly graded offset."""
    min_offset = 0.0
    max_offset = 3.0
    radius = 0.5 * (np.max(x) - np.min(x))
    return (max_offset - min_offset) * (x**2 + y**2) / radius**2 + min_offset


geometry = Tpms(
    surface_function=gyroid,
    offset=linear_graded_offset,
    repeat_cell=(5, 2, 1),
    resolution=30,
)
sheet = geometry.sheet

plotter = pv.Plotter()
plotter.add_mesh(geometry.sheet, color="w")
plotter.view_xy()
plotter.enable_parallel_projection()
plotter.show_axes()
# plotter.show()
