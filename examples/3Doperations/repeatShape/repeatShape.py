import cadquery as cq
from microgen import Rve, repeatShape, Phase

rve = Rve(dim_x=1, dim_y=1, dim_z=1)

unit_geom = cq.importers.importStep(
    "../../Lattices/octetTruss/octettruss.step"
)

unit_geom = cq.Shape(unit_geom.val().wrapped)

shape = repeatShape(unit_geom, rve, grid=(5, 5, 5))

cq.exporters.export(shape, "repeated_shape.stl")
