from pathlib import Path

import cadquery as cq

from microgen import Rve, repeatShape

rve = Rve(dim=1)

step_file = str(Path(__file__).parent / "octettruss.step")
unit_geom = cq.importers.importStep(step_file)

unit_geom = cq.Shape(unit_geom.val().wrapped)

shape = repeatShape(unit_geom, rve, grid=(5, 5, 5))

stl_file = str(Path(__file__).parent / "repeated_shape.stl")
cq.exporters.export(shape, stl_file)
