# from microgen import Rve, Tpms, Sphere, repeatGeometry, tpms, Phase
# import cadquery as cq

# rve = Rve(1, 1, 1)
# geometry = Tpms(
#     center=(0.5, 0.5, 0.5),
#     surface_function=tpms.gyroid,
#     type_part="sheet",
#     thickness=0.3,
#     path_data="data",
#     cell_size=1.,
#     repeat_cell=3
# )
# shape = geometry.generate()

# sphere = Sphere(center=(1, 1, 1), radius=1.5)
# sphere_shape = sphere.generate()
# print("ok")

# result = shape.intersect(sphere_shape)
# print("ok")
# cq.exporters.export(result, "tpms_sphere.stl")

from microgen import Rve, Tpms, Sphere, repeatGeometry, tpms, Phase
import cadquery as cq

rve = Rve(1, 1, 1)
geometry = Tpms(
    center=(0.5, 0.5, 0.5),
    surface_function=tpms.gyroid,
    type_part="sheet",
    thickness=0.2,
    path_data="data",
)
shape = geometry.generate(sizeMesh=0.03, minFacetAngle=20.0, maxRadius=0.03)

print("repeat geometry")
phase = repeatGeometry(Phase(shape=shape), rve, grid=(3, 3, 3))

sphere = Sphere(center=(1, 1, 1), radius=1.5)
sphere_shape = sphere.generate()

result = phase.shape.intersect(sphere_shape)
cq.exporters.export(result, "tpms_sphere.stl")