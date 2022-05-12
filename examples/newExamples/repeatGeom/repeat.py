import cadquery as cq
from microgen import Rve, fuseParts, BasicGeometry, repeatGeometry

# Size of the mesh
size_mesh = 0.03
a = 1.0
b = 1.0
c = 1.0

rve = Rve(dim_x=a, dim_y=b, dim_z=c, size_mesh=size_mesh)

# unit_geom = cq.importers.importStep('../../octetTruss/compound.step')
unit_geom = cq.importers.importStep('../../triplyPeriodicMinimalSurfaces/sheet.step')
unit_geom = cq.Shape(unit_geom.val().wrapped)

final_geom = repeatGeometry(unit_geom, rve, grid={"x": 5, 
                                                  "y": 5, 
                                                  "z": 5})

cq.exporters.export(final_geom, 'compound.stl')