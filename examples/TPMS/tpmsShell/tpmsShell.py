from microgen import Tpms, tpms
import cadquery as cq

shell = (
    cq.Workplane("front")
    .box(3, 3, 3)
    .faces("+Z or -X or +X")
    .shell(0.1)
)
shell = shell.val()

geometry = Tpms(
    surface_function=tpms.gyroid,
    type_part="sheet",
    thickness=0.2,
    repeat_cell=3,
)
shape = geometry.generate(nSample = 30, smoothing=0)

compound = cq.Compound.makeCompound([shell, shape])
cq.exporters.export(compound, "tpms_shell.stl")
