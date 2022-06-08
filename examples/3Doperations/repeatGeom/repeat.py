import cadquery as cq
from microgen import Rve, repeatGeometry, Phase

rve = Rve(dim_x=1, dim_y=1, dim_z=1)

unit_geom = cq.importers.importStep('../../3Doperations/BooleanOperations/octetTruss/octettruss.step')
# unit_geom = cq.importers.importStep('../../TPMS/Gyroid/sheet.step')

unit_geom = cq.Shape(unit_geom.val().wrapped)

final_geom = repeatGeometry(Phase(shape=unit_geom), rve, grid={"x": 5, "y": 5, "z": 5})

cq.exporters.export(final_geom, 'repeated_geometry.stl')
