from pathlib import Path

import pyvista as pv

from microgen import Tpms, surface_functions

geometry = Tpms(
    surface_function=surface_functions.gyroid,
    offset=0.3,
    repeat_cell=3,
    resolution=30,
)
shape = geometry.generateVtk(type_part="sheet")
shape.flip_normals()

sphere = pv.Sphere(radius=1.45)

result = shape.boolean_intersection(sphere)
vtk_file = Path(__file__).parent / "tpmsSphere.vtk"
result.save(vtk_file)
# result.plot(color='w')
