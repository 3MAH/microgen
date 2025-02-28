from functools import partial

import numpy as np
from scipy.spatial.transform import Rotation

from microgen import CylindricalTpms
from microgen.shape.surface_functions import gyroid


def rotated_gyroid(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    grid_rotation: Rotation,
) -> np.ndarray:
    x_rot, y_rot, z_rot = grid_rotation.apply(np.stack([x, y, z], axis=-1)).T
    return gyroid(x=x_rot, y=y_rot, z=z_rot)


grid_rotation = Rotation.from_euler("z", 45, degrees=True)

geometry = CylindricalTpms(
    surface_function=partial(rotated_gyroid, grid_rotation=grid_rotation),
    offset=0.5,
    cell_size=(1, 1, 1),
    repeat_cell=(1, 0, 1),
    radius=1,
    resolution=20,
)
sheet = geometry.sheet

# import pyvista as pv
# plotter = pv.Plotter()
# plotter.add_mesh(sheet, color="w")
# plotter.view_xy()
# plotter.enable_parallel_projection()
# plotter.show_axes()
# plotter.show()
