from pathlib import Path

import numpy as np

from microgen.cad import make_compound
from microgen.shape import new_geometry

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


parts = []

n_col = 3
n_row = np.ceil(len(shapes) / n_col)
i = 0
for shape, param_geom in shapes.items():
    i_x = i % n_col
    i_y = i // n_col
    elem = new_geometry(
        shape=shape,
        center=(1.2 * (i_x - 0.5 * (n_col - 1)), -1.2 * (i_y - 0.5 * (n_row - 1)), 0),
        orientation=(90, 90, 90),
        param_geom=param_geom,
    )
    parts.append(elem.generate())
    i = i + 1

compound = make_compound(parts)
stl_file = str(Path(__file__).parent / "shapes.stl")
compound.export_stl(stl_file)
