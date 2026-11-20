"""
Spherical gyroid
================

Intersection of a gyroid sheet with a sphere, computed with a boolean
operation on the surface meshes.
"""

import pyvista as pv

from microgen import Tpms, surface_functions

geometry = Tpms(
    surface_function=surface_functions.gyroid,
    offset=0.3,
    repeat_cell=3,
    resolution=30,
)
shape = geometry.generate_surface_mesh(type_part="sheet")
shape = shape.flip_faces()

sphere = pv.Sphere(radius=1.45)

result = shape.boolean_intersection(sphere)
result.save("tpms_sphere.vtk")

# %%
result.plot(smooth_shading=True, color="#ef7d2d")
