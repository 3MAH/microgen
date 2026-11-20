"""
Basic shapes
============

Create the elementary shapes of microgen with :func:`microgen.shape.new_geometry`
and export them as a single CAD compound.

Requires the CAD backend: ``pip install microgen[cad]``.
"""

import numpy as np
import pyvista as pv

from microgen.cad import make_compound
from microgen.shape import new_geometry

# %%
# Shape parameters, keyed by shape name.
shapes = {
    "box": {"dim": (0.8, 0.8, 0.8)},
    "sphere": {"radius": 0.5},
    "capsule": {"height": 0.5, "radius": 0.3},
    "cylinder": {"height": 0.5, "radius": 0.5},
    "ellipsoid": {"radii": (0.5, 0.25, 0.3)},
    "extrudedpolygon": {
        "list_corners": [
            (0.5, 0),
            (0.25, 0.44),
            (-0.25, 0.44),
            (-0.5, 0),
            (-0.25, -0.44),
            (0.25, -0.44),
            (0.5, 0),
        ],
        "height": 0.5,
    },
}

# %%
# Place the shapes on a grid and build their CAD representation.
parts = []
n_col = 3
n_row = np.ceil(len(shapes) / n_col)
for i, (shape, param_geom) in enumerate(shapes.items()):
    i_x = i % n_col
    i_y = i // n_col
    elem = new_geometry(
        shape=shape,
        center=(1.2 * (i_x - 0.5 * (n_col - 1)), -1.2 * (i_y - 0.5 * (n_row - 1)), 0),
        orientation=(90, 90, 90),
        param_geom=param_geom,
    )
    parts.append(elem.generate_cad())

compound = make_compound(parts)
stl_file = "shapes.stl"
compound.export_stl(stl_file)

# %%
# Display the exported compound.
pv.read(stl_file).plot(smooth_shading=True, split_sharp_edges=True, color="#ef7d2d")
