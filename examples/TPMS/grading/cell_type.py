"""Example of morphing between multiple TPMS surface types.

The morphing between TPMS is managed by the TPMS_TYPES dictionary that maps
control points to TPMS functions. The weight function is used to determine the
contribution of each TPMS function to the final surface.

**Note:** The weight function is an exponential function that transitions
between control points.
"""

from __future__ import annotations

from typing import Callable, Tuple

import numpy as np
import numpy.typing as npt
import pyvista as pv

from microgen import Tpms
from microgen.shape.surface_functions import fischerKochS, gyroid, schwarzD, schwarzP

ControlPoint = Tuple[float, float, float]

REPEAT = (4, 4, 2)
TPMS_TYPES: dict[ControlPoint, Callable] = {
    (-1, 0, 0): schwarzP,
    (1, 0, 0): gyroid,
    (0, 1, 0): fischerKochS,
    (0, -1, 0): schwarzD,
}


def exp_func(
    x: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
    z: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Transition between points using exponential function."""
    k = 0.5
    norm = x**2 + y**2 + z**2
    return 1 + np.exp(k * norm)


def weight(
    x: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
    z: npt.NDArray[np.float64],
    control_points: list[ControlPoint],
    index: int,
) -> npt.NDArray[np.float64]:
    """Weight function for morphing between points."""
    denom = np.zeros_like(x)
    for p in control_points:
        denom += exp_func(x + p[0], y + p[1], z + p[2])
    point = control_points[index]
    return exp_func(x + point[0], y + point[1], z + point[2]) / denom


def multi_morph(
    tpms_types: dict[ControlPoint, Callable],
    x: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
    z: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Morphing between multiple functions."""
    control_points = list(tpms_types.keys())
    surfaces = list(tpms_types.values())
    result = np.zeros_like(x)
    for index, surface_function in enumerate(surfaces):
        weight_func = weight(x, y, z, control_points, index)
        result += weight_func * surface_function(x, y, z)
    return result


def trigraded(
    x: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
    z: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Morph between three TPMS functions."""
    return multi_morph(
        tpms_types=TPMS_TYPES,
        x=x,
        y=y,
        z=z,
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
# plotter.show()
