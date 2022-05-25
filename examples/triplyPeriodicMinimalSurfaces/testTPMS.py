import cadquery as cq
from microgen import Rve, BasicGeometry, gyroid

a = 1.0
b = 1.0
c = 1.0

periodicity = 0
rve = Rve(a, b, c)

elem = BasicGeometry(number=0, shape='tpms',
                     xc=0.5, yc=0.5, zc=0.5,
                     psi=0., theta=0., phi=0.,
                     param_geom={"surface_function": gyroid,
                                 "type_part": 'sheet',
                                 "thickness": 0.1},
                     path_data='data')
elem.geometry.createSurfaces(rve=rve,
                             sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03,
                             path_data='data')
part = elem.generate(rve=rve)

cq.exporters.export(part, 'gyroid.step')
