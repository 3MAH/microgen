from microgen import Tpms
from microgen.shape.surface_functions import gyroid

geometry = Tpms(
    surface_function=gyroid,
    density=0.30,
    resolution=30,
)
shape = geometry.generateVtk(type_part="sheet").extract_surface()
shape.save("gyroid.stl")
