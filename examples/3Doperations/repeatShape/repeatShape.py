import cadquery as cq

from microgen import Rve, repeatShape

rve = Rve(dim=1)

unit_geom = cq.importers.importStep("octettruss.step")

unit_geom = cq.Shape(unit_geom.val().wrapped)

shape = repeatShape(unit_geom, rve, grid=(5, 5, 5))

cq.exporters.export(shape, "repeated_shape.stl")
