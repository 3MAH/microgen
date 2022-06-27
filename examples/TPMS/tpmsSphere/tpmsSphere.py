from microgen import Tpms, Sphere, tpms
import cadquery as cq

geometry = Tpms(
    surface_function=tpms.gyroid,
    type_part="sheet",
    thickness=0.4,
    repeat_cell=3
)
shape = geometry.generate()
cq.exporters.export(shape, "gyroid.stl")
print('gyroid')

sphere = Sphere(radius=1.45)
sphere_shape = sphere.generate()
cq.exporters.export(sphere_shape, "sphere.stl")
print('sphere')

result = shape.intersect(sphere_shape)
cq.exporters.export(result, "tpms_sphere.stl")
print('intersect')
