"""Example of morphing between multiple TPMS surface types."""

from typing import Callable, List, Tuple

import numpy as np
import pyvista as pv

from microgen import Tpms
from microgen.shape.surface_functions import fischerKochS, gyroid, schwarzD, schwarzP

REPEAT = (4, 4, 2)
POINTS = [(-0.5, 0, 0), (0.5, 0, 0), (0, 0.5, 0), (0, -0.5, 0)]
SURFACES = [schwarzP, gyroid, fischerKochS, schwarzD]


def exp_func(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Transition between points using exponential function."""
    k = 0.5
    norm = x**2 + y**2 + z**2
    return 1 + np.exp(k * norm)


def weight(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    points: List[Tuple[float, float, float]],
    index: int,
) -> np.ndarray:
    """Weight function for morphing between points."""
    denom = np.zeros_like(x)
    for p in points:
        denom += exp_func(x - p[0], y - p[1], z - p[2])
    point = points[index]
    return exp_func(x - point[0], y - point[1], z - point[2]) / denom


def multi_morph(
    phi: List[Callable],
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    points: List[Tuple[float, float, float]],
) -> np.ndarray:
    """Morphing between multiple functions."""
    result = np.zeros_like(x)
    for index, surface_function in enumerate(phi):
        weight_func = weight(x, y, z, points, index)
        result += weight_func * surface_function(x, y, z)
    return result


def trigraded(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Morph between three TPMS functions."""
    return multi_morph(
        phi=SURFACES,
        x=x,
        y=y,
        z=z,
        points=POINTS,
    )


geometry = Tpms(
    surface_function=trigraded,
    offset=0.3,
    repeat_cell=REPEAT,
    resolution=50,
)
sheet = geometry.sheet

plotter = pv.Plotter()
plotter.add_mesh(geometry.sheet, color="w")
plotter.view_xy()
plotter.enable_parallel_projection()
plotter.show_axes()
plotter.show()
