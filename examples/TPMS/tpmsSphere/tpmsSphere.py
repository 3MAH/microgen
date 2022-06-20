from microgen import Tpms, Sphere, tpms
import cadquery as cq

# unit_geom = Tpms(
#     surface_function=tpms.gyroid,
#     type_part="sheet",
#     thickness=0.2,
#     path_data="data",
# )
# cq.exporters.export(unit_geom.generate(), "unit.stl")

geometry = Tpms(
    surface_function=tpms.gyroid,
    type_part="sheet",
    thickness=0.2,
    cell_size=1/2,
    repeat_cell=2,
    path_data="data",
)
shape = geometry.generate()

sphere = Sphere()
sphere_shape = sphere.generate()

result = shape.intersect(sphere_shape)
cq.exporters.export(result, "tpms_sphere.stl")
