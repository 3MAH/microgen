from microgen import Tpms, surface_functions
import cadquery as cq

shell = (
    cq.Workplane("front")
    .box(3, 3, 3)
    .faces("+Z or -X or +X")
    .shell(0.1)
)
shell = shell.val()

geometry = Tpms(
    surface_function=surface_functions.gyroid,
    offset=0.5,
    repeat_cell=3,
    resolution=20,
)
shape = geometry.generate(type_part="sheet", smoothing=0)

compound = cq.Compound.makeCompound([shell, shape])
cq.exporters.export(compound, "tpms_shell.stl")
