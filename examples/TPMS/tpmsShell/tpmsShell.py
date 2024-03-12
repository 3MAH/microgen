from pathlib import Path

import cadquery as cq

from microgen import Tpms, surface_functions

shell = cq.Workplane("front").box(3, 3, 3).faces("+Z or -X or +X").shell(0.1)
shell = shell.val()

geometry = Tpms(
    surface_function=surface_functions.gyroid,
    offset=0.5,
    repeat_cell=3,
    resolution=15,
)
shape = geometry.generate(type_part="sheet", smoothing=0)

compound = cq.Compound.makeCompound([shell, shape])
stl_file = str(Path(__file__).parent / "tpms_shell.stl")
cq.exporters.export(compound, stl_file)
