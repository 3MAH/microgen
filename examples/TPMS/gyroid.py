"""
Gyroid
======

A gyroid sheet of relative density 0.3, exported as an STL surface mesh.
"""

from microgen import Tpms
from microgen.shape.surface_functions import gyroid

geometry = Tpms(
    surface_function=gyroid,
    density=0.30,
    resolution=30,
)
shape = geometry.generate_surface_mesh(type_part="sheet")
shape.save("gyroid.stl")

# %%
shape.plot(smooth_shading=True, color="#ef7d2d")
