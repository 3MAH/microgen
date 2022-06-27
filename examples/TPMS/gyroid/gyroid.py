from microgen import Tpms, tpms
import cadquery as cq

geometry = Tpms(
    surface_function=tpms.gyroid,
    type_part="sheet",
    thickness=0.2
)
shape = geometry.generate()

cq.exporters.export(shape, "gyroid.stl")
