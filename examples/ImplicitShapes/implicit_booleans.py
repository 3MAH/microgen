"""
Boolean operations on implicit fields
=====================================

Since microgen 2.0, every shape carries an implicit field ``f(x, y, z)``,
negative inside. Booleans are then operations on numbers, not on
boundary representations:

- union ``a | b`` is ``min(fa, fb)``,
- intersection ``a & b`` is ``max(fa, fb)``,
- difference ``a - b`` is ``max(fa, -fb)``,
- complement ``~a`` is ``-fa``.

They cost nothing until the shape is meshed, never fail on tangent or
coincident faces, apply to TPMS, spinodoids, lattices and primitives
alike, and need no CAD kernel. Smooth variants add fillets of radius ``k``
at the junctions, which are hard to get with a BREP kernel.
"""

# sphinx_gallery_thumbnail_number = 1
import time

import pyvista as pv

from microgen import Box, Cylinder, Sphere, Tpms, implicit_ops, surface_functions

# %%
# Constructive solid geometry
# ---------------------------
# The classic CSG example: a box intersected with a sphere, minus three
# cylinders along the axes.
box = Box(dim=(1.0, 1.0, 1.0))
sphere = Sphere(radius=0.65)
cylinders = (
    Cylinder(height=1.5, radius=0.3)
    | Cylinder(height=1.5, radius=0.3, orientation=(90.0, 0.0, 0.0))
    | Cylinder(height=1.5, radius=0.3, orientation=(0.0, 90.0, 90.0))
)

start = time.perf_counter()
csg = (box & sphere) - cylinders
print(f"Composing the CSG tree: {1e3 * (time.perf_counter() - start):.2f} ms")

bounds = (-0.6, 0.6, -0.6, 0.6, -0.6, 0.6)
start = time.perf_counter()
csg_mesh = csg.generate_surface_mesh(bounds=bounds, resolution=120)
print(f"Meshing it: {time.perf_counter() - start:.2f} s, {csg_mesh.n_cells} triangles")

# %%
# Smooth booleans
# ---------------
# ``smooth_union``, ``smooth_intersection`` and ``smooth_difference`` blend
# the two fields over a distance ``k``: the edges of the CSG part become
# fillets.
smooth = box.smooth_intersection(sphere, k=0.1).smooth_difference(cylinders, k=0.1)
smooth_mesh = smooth.generate_surface_mesh(bounds=bounds, resolution=120)
smooth_mesh.plot(color="#ef7d2d", smooth_shading=True)

# %%
# Mixing primitives and TPMS
# --------------------------
# The same operators combine any implicit shapes. Here, a gyroid is
# clipped by a sphere and joined to a solid core with a smooth union.
# :func:`microgen.implicit_ops.shell` turns a solid into a hollow wall of
# given thickness.
gyroid = Tpms(
    surface_function=surface_functions.gyroid, offset=0.6, cell_size=0.5, repeat_cell=3
)
gyroid_ball = (gyroid & Sphere(radius=0.7)).smooth_union(Sphere(radius=0.25), k=0.1)
gyroid_ball_mesh = gyroid_ball.generate_surface_mesh(
    bounds=(-0.75, 0.75, -0.75, 0.75, -0.75, 0.75), resolution=150
)

octant = Box(dim=(1.0, 1.0, 1.0), center=(0.5, 0.5, 0.5))
shelled = implicit_ops.shell(Sphere(radius=0.6), thickness=0.08) - octant
shelled_mesh = shelled.generate_surface_mesh(bounds=bounds, resolution=120)

# %%
# Hard CSG, smooth CSG, a gyroid ball with a solid core, and a hollow sphere
# cut open by removing one octant.
pl = pv.Plotter(shape=(2, 2), window_size=(1000, 1000))
for i, (title, mesh) in enumerate(
    [
        ("(box & sphere) - cylinders", csg_mesh),
        ("smooth booleans, k = 0.1", smooth_mesh),
        ("gyroid & sphere, smooth union with a core", gyroid_ball_mesh),
        ("shell(sphere) - octant", shelled_mesh),
    ]
):
    pl.subplot(i // 2, i % 2)
    pl.add_text(title, font_size=10)
    pl.add_mesh(mesh, color="#ef7d2d", smooth_shading=True)
pl.link_views()
pl.show()
