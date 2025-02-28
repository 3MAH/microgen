import numpy as np
from functools import partial
from scipy.spatial.transform import Rotation as R
from microgen import CylindricalTpms
from microgen.shape.surface_functions import gyroid


def rotated_gyroid(x, y, z, grid_rotation_matrix):
    coords_rotated = np.dot(
        np.stack([x, y, z], axis=-1), (grid_rotation_matrix.as_matrix()).T
    )
    x_rot, y_rot, z_rot = coords_rotated.T

    return gyroid(x=x_rot, y=y_rot, z=z_rot)


grid_rotation_matrix = R.from_euler("zxz", [45.0, 0, 0], degrees=True)

geometry = CylindricalTpms(
    surface_function=partial(rotated_gyroid, grid_rotation_matrix=grid_rotation_matrix),
    offset=0.5,
    cell_size=(1, 1, 1),
    repeat_cell=(1, 0, 1),
    radius=1,
    resolution=20,
)
sheet = geometry.sheet

# plotter = pv.Plotter()
# plotter.add_mesh(sheet, color="w")
# plotter.view_xy()
# plotter.enable_parallel_projection()
# plotter.show_axes()
# plotter.show()
