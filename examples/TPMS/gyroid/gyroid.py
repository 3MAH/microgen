from microgen import Tpms
from microgen.shape.surface_functions import gyroid

geometry = Tpms(surface_function=gyroid).with_density(0.30).with_resolution(30)
shape = geometry.generate_surface_mesh(type_part="sheet")
shape.save("gyroid.stl")
