from microgen import Tpms
from microgen.shape.surface_functions import gyroid

geometry = Tpms(
    surface_function=gyroid,
    offset=0.3,
    resolution=100,
)
shape = geometry.generateVtk(type_part="sheet")
shape.save("gyroid.stl")
