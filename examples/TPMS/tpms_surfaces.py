"""
Available TPMS
==============

The surface functions available in :mod:`microgen.shape.surface_functions`:
Neovius, gyroid, Schwarz D, Schwarz P, honeycomb, Schoen IWP, Schoen FRD,
Fischer-Koch S and PMY.
"""

import numpy as np
import pyvista as pv

from microgen import Tpms
from microgen.shape import surface_functions

surfaces = [
    surface_functions.neovius,
    surface_functions.gyroid,
    surface_functions.schwarz_d,
    surface_functions.schwarz_p,
    surface_functions.honeycomb,
    surface_functions.schoen_iwp,
    surface_functions.schoen_frd,
    surface_functions.fischer_koch_s,
    surface_functions.pmy,
]

meshes = []
n_col = 3
n_row = np.ceil(len(surfaces) / n_col)
for i, surface in enumerate(surfaces):
    i_x = i % n_col
    i_y = i // n_col
    elem = Tpms(
        surface_function=surface,
        offset=0.3,
        resolution=50,
    )
    mesh = elem.sheet
    mesh.translate(
        [1.2 * (i_x - 0.5 * (n_col - 1)), -1.2 * (i_y - 0.5 * (n_row - 1)), 0],
        inplace=True,
    )
    meshes.append(mesh)

surfaces_block = pv.MultiBlock(meshes)
surfaces_block.save("tpms_surfaces.vtm")

# %%
surfaces_block.plot(smooth_shading=True, cpos="xy", color="#ef7d2d")
