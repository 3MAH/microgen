"""
Fischer-Koch S surface
======================

The zero level set of the Fischer-Koch S function on a 5 x 5 x 5 array of cells.
"""

from microgen import Tpms
from microgen.shape.surface_functions import fischer_koch_s

geometry = Tpms(surface_function=fischer_koch_s, repeat_cell=5, offset=0)
mesh = geometry.generate_surface_mesh(type_part="surface")
mesh.save("fischer_koch_s.vtk")

# %%
mesh.plot(smooth_shading=True, color="#ef7d2d")
