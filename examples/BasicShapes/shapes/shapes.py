from pathlib import Path

import cadquery as cq
import numpy as np

from microgen.shape import newGeometry

shapes = {
    "box": {"dim_x": 0.8, "dim_y": 0.8, "dim_z": 0.8},
    "sphere": {"radius": 0.5},
    "capsule": {"height": 0.5, "radius": 0.3},
    "cylinder": {"height": 0.5, "radius": 0.5},
    "ellipsoid": {"a_x": 0.5, "a_y": 0.25, "a_z": 0.3},
    "extrudedpolygon": {
        "listCorners": [
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


assembly = cq.Assembly()

n_col = 3
n_row = np.ceil(len(shapes) / n_col)
i = 0
for shape, param_geom in shapes.items():
    i_x = i % n_col
    i_y = i // n_col
    elem = newGeometry(
        shape=shape,
        center=(1.2 * (i_x - 0.5 * (n_col - 1)), -1.2 * (i_y - 0.5 * (n_row - 1)), 0),
        orientation=(90, 90, 90),
        param_geom=param_geom,
    )
    assembly.add(elem.generate())
    i = i + 1

compound = assembly.toCompound()
stl_file = str(Path(__file__).parent / "shapes.stl")
cq.exporters.export(compound, stl_file)
