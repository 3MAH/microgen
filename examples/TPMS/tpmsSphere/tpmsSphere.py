from pathlib import Path

import pyvista as pv

from microgen import Tpms, surface_functions

geometry = (
    Tpms(surface_function=surface_functions.gyroid)
    .with_offset(0.3)
    .with_repeat_cell(3)
    .with_resolution(30)
)
shape = geometry.generate_surface_mesh(type_part="sheet")
shape = shape.flip_faces()

sphere = pv.Sphere(radius=1.45)

result = shape.boolean_intersection(sphere)
vtk_file = Path(__file__).parent / "tpmsSphere.vtk"
result.save(vtk_file)
# result.plot(color='w')
