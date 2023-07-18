from microgen import Tpms, surface_functions
import pyvista as pv

geometry = Tpms(
    surface_function=surface_functions.gyroid,
    offset=0.3,
    repeat_cell=3,
    resolution=30,
)
shape = geometry.generateVtk(type_part="sheet")
shape.triangulate(inplace=True)
shape.flip_normals()

sphere = pv.Sphere(radius=1.45)

result = shape.boolean_intersection(sphere)
result.plot(color='w')
