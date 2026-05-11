from pathlib import Path

from microgen import Rve, repeat_shape
from microgen.cad import import_step

rve = Rve(dim=1)

step_file = str(Path(__file__).parent / "octettruss.step")
unit_geom = import_step(step_file)

shape = repeat_shape(unit_geom, rve, grid=(5, 5, 5))

stl_file = str(Path(__file__).parent / "repeated_shape.stl")
shape.export_stl(stl_file)
