import cadquery as cq
from microgen import Rve, fuseParts, BasicGeometry, repeatGeometry

a = 1.0
b = 1.0
c = 1.0

rve = Rve(dim_x=a, dim_y=b, dim_z=c)

unit_geom = cq.importers.importStep('../../octetTruss/octettruss.step')
# unit_geom = cq.importers.importStep('../../triplyPeriodicMinimalSurfaces/sheet.step')
unit_geom = cq.Shape(unit_geom.val().wrapped)

final_geom = repeatGeometry(unit_geom, rve, grid={"x": 5, 
                                                  "y": 5, 
                                                  "z": 5})

cq.exporters.export(final_geom, 'repeated_geometry.stl')