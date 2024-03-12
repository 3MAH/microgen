from pathlib import Path

import numpy as np
import pyvista as pv

from microgen import Tpms
from microgen.shape import surface_functions

surfaces = [
    surface_functions.neovius,
    surface_functions.gyroid,
    surface_functions.schwarzD,
    surface_functions.schwarzP,
    surface_functions.honeycomb,
    surface_functions.schoenIWP,
    surface_functions.schoenFRD,
    surface_functions.fischerKochS,
    surface_functions.pmy,
]

meshes = []

i = 0
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

    # center = (1.2 * (i_x - 0.5 * (n_col - 1)), -1.2 * (i_y - 0.5 * (n_row - 1)), 0)
    mesh = elem.sheet
    mesh.translate(
        [1.2 * (i_x - 0.5 * (n_col - 1)), -1.2 * (i_y - 0.5 * (n_row - 1)), 0],
        inplace=True,
    )
    meshes.append(mesh)

vtm_file = Path(__file__).parent / "surfaces.vtm"
pv.MultiBlock(meshes).save(vtm_file)
# pv.MultiBlock(meshes).plot(color="w")
