from microgen import Rve, Tpms, Sphere, repeatGeometry, tpms, Phase
import cadquery as cq

rve = Rve(1, 1, 1)
geometry = Tpms(center=(0.5, 0.5, 0.5),
                surface_function=tpms.gyroid,
                type_part='sheet',
                thickness=0.2,
                path_data='data')
shape = geometry.generate(sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03)

print("repeat geometry")
compound = repeatGeometry(Phase(shape=shape), rve, {'x': 3, 'y': 3, 'z': 3})

sphere = Sphere(center=(1, 1, 1), radius=1.5)
sphere_shape = sphere.generate()

result = compound.intersect(sphere_shape)
cq.exporters.export(result, 'tpms_sphere.stl')
