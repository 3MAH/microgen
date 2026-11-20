"""
Gyroid in a shell
=================

A gyroid sheet inside a hollow box, as one implicit shape: the hollow box is
a difference of two boxes, joined to the gyroid sheet with a union. One
octant is removed to show the inside.

The result is meshed by marching cubes, without the CAD backend. With
``pip install microgen[cad]``, ``Phase.from_shape(part).cad`` gives the
OCCT shape for a STEP export.
"""

from microgen import Box, Tpms, implicit_ops, surface_functions

# %%
# Outer 3 x 3 x 3 box, hollowed by subtracting an inner 2.8 x 2.8 x 2.8 cavity.
hollow_box = implicit_ops.box((3.0, 3.0, 3.0)) - implicit_ops.box((2.8, 2.8, 2.8))

# %%
# 27 gyroid cells filling the box, as a sheet 0.1 thick, clipped by the box.
gyroid = Tpms(surface_function=surface_functions.gyroid, offset=0.3, repeat_cell=3)
sheet = gyroid.as_sheet(thickness=0.1) & implicit_ops.box((2.9, 2.9, 2.9))
part = hollow_box | sheet

# %%
# Remove the octant facing the camera and mesh the field.
octant = Box(dim=(1.6, 1.6, 1.6), center=(0.8, 0.8, 0.8))
cut_part = part - octant
mesh = cut_part.generate_surface_mesh(
    bounds=(-1.55, 1.55, -1.55, 1.55, -1.55, 1.55), resolution=150
)
mesh.save("tpms_shell.stl")
print(f"{mesh.n_cells} triangles")

mesh.plot(color="#ef7d2d", smooth_shading=True)
