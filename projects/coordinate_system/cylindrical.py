import numpy as np
import pyvista as pv

from microgen import pv_tpms


def swapped_gyroid(x, y, z):
    return np.cos(x) * np.sin(y) + np.cos(y) * np.sin(z) + np.cos(z) * np.sin(x)

geometry = pv_tpms.Tpms(surface_function=swapped_gyroid,
                        offset=0.5,
                        cell_size=(1, 1, 1),
                        repeat_cell=(1, 6, 1),
                        coordinate_system="cylindrical",
                        cylinder_radius=1)

pl = pv.Plotter()
pl.add_mesh(geometry.sheet, style="surface", color="white")
pl.parallel_projection = True
pl.view_xy()
pl.show_axes()
pl.show_grid()
pl.show()
