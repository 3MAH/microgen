import microgen
import cadquery as cq

rve = microgen.Rve(1, 1, 1)
phase = microgen.BasicGeometry(number=0, shape='tpms',
                               xc=0.5, yc=0.5, zc=0.5,
                               psi=0., theta=0., phi=0.,
                               param_geom={"surface_function": microgen.shape.tpms.schwarzP,
                                           "type_part": 'sheet',
                                           "thickness": 0.2},
                               path_data='data')
phase.geometry.createSurfaces(rve=rve,
                              sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03,
                              path_data='data')
gyroid = phase.generate(rve=rve)

print("repeat geometry")
gyroid = microgen.repeatGeometry(gyroid, rve, {'x': 3, 'y': 3, 'z': 3})
# cq.exporters.export(gyroid, 'gyroid.stl')

sphere = microgen.BasicGeometry(number=0, shape='sphere',
                                xc=1, yc=1, zc=1,
                                psi=0., theta=0., phi=0.,
                                param_geom={"radius": 1.5})
sphere = sphere.generate()
# cq.exporters.export(sphere, 'sphere.stl')

result = gyroid.intersect(sphere)
cq.exporters.export(result, 'tpms_sphere.stl')
