from typing import List

import numpy as np
import pyvista as pv

from microgen import pv_tpms


def swapped_gyroid(x, y, z):
    return np.cos(x) * np.sin(y) + np.cos(y) * np.sin(z) + np.cos(z) * np.sin(x)


geometry = pv_tpms.Tpms(surface_function=swapped_gyroid,
                        offset=1.,)

print(f"relative density = {geometry.sheet.volume / geometry.grid.volume:.2%}")

pl = pv.Plotter()
pl.add_mesh(geometry.sheet, color="w")
pl.add_mesh(geometry.upper_skeletal, color="b")
pl.add_mesh(geometry.lower_skeletal, color="r")
pl.show()
