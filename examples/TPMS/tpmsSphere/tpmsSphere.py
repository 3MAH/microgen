from microgen import Tpms, surface_functions
import pyvista as pv

geometry = Tpms(
    surface_function=surface_functions.gyroid,
    offset=0.3,
    repeat_cell=3,
    resolution=30,
)
shape = geometry.generateVtk(type_part="sheet").extract_surface()
shape.triangulate(inplace=True)
shape.flip_normals()

sphere = pv.Sphere(radius=1.45)

result = shape.boolean_intersection(sphere)
result.save("tpmsSphere.vtk")
# result.plot(color='w')
