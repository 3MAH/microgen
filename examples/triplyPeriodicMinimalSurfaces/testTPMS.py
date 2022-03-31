import cadquery as cq
from microgen import Rve, BasicGeometry

# Size of the mesh
size_mesh = 0.03
a = 1.0
b = 1.0
c = 1.0

periodicity = 0
revel = Rve(a, b, c, size_mesh)
phases = []

elem = BasicGeometry(number=0, shape='tpms',
                     xc=0.5, yc=0.5, zc=0.5,
                     psi=0., theta=0., phi=0.,
                     param_geom=['gyroid', 'sheet', 0.1], path_data='data')
elem.geometry.createSurfaces(rve=revel,
                             sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03,
                             path_data='data')
sheet = elem.generate(rve=revel)

cq.exporters.export(sheet, 'sheet.step')
