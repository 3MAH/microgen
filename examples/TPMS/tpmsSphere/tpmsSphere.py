from microgen import Tpms, Sphere, tpms
import cadquery as cq

geometry = Tpms(
    surface_function=tpms.gyroid,
    type_part="sheet",
    thickness=0.2,
    repeat_cell=3,
    path_data="data",
)
shape = geometry.generate(sizeMesh=0.05, minFacetAngle=20., maxRadius=0.05)
cq.exporters.export(shape, "gyroid.stl")
print('gyroid')

sphere = Sphere(radius=1.2)
sphere_shape = sphere.generate()
cq.exporters.export(sphere_shape, "sphere.stl")
print('sphere')

result = shape.intersect(sphere_shape)
cq.exporters.export(result, "tpms_sphere.stl")
print('intersect')
