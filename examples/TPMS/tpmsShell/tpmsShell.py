from pathlib import Path

from microgen import Tpms, surface_functions
from microgen.cad import make_box, make_compound

# Outer 3x3x3 box, hollowed by subtracting an inner 2.8x2.8x2.8 cavity —
# replaces the CadQuery ``Workplane.shell`` API.  The face-selective version
# ("open on +Z, -X, +X") could be reproduced with ``BRepOffsetAPI_MakeThickSolid``
# via OCP, but a plain symmetric shell is enough for this demo.
outer = make_box((3.0, 3.0, 3.0), (0.0, 0.0, 0.0))
inner = make_box((2.8, 2.8, 2.8), (0.0, 0.0, 0.0))
shell = outer.cut(inner)

geometry = (
    Tpms(surface_function=surface_functions.gyroid)
    .with_offset(0.5)
    .with_repeat_cell(3)
    .with_resolution(15)
)
shape = geometry.generate_cad(type_part="sheet", smoothing=0)

compound = make_compound([shell, shape])
stl_file = str(Path(__file__).parent / "tpms_shell.stl")
compound.export_stl(stl_file)
