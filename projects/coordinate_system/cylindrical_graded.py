import numpy as np
import pyvista as pv

from microgen import pv_tpms


def swapped_gyroid(x, y, z):
    return np.cos(x) * np.sin(y) + np.cos(y) * np.sin(z) + np.cos(z) * np.sin(x)

def offset(x, y, z):
    a = (hmax - hmin) / (np.max(x) - np.min(x))
    b = hmax - a * np.max(x)

    return a * x + b

geometry = pv_tpms.Tpms(surface_function=swapped_gyroid,
                        offset=offset,
                        cell_size=(1, 1, 1),
                        repeat_cell=(10, 1000, 1),
                        coordinate_system="cylindrical",
                        cylinder_radius=10)


hmin = 0.
hmax = 3.

pl = pv.Plotter()
pl.add_mesh(geometry.sheet, style="surface", color="white")
pl.parallel_projection = True
pl.view_xy()
pl.show_axes()
pl.show_grid()
pl.show()
