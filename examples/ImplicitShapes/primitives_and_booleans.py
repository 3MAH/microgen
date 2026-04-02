"""Implicit shapes: primitives, boolean operations, and smooth blending.

Demonstrates the F-rep (Function Representation) framework in microgen
for creating implicit shapes and combining them with hard and smooth
boolean operations.
"""

import pyvista as pv

from microgen.shape.implicit_shape import (
    batch_smooth_union,
    implicit_box,
    implicit_sphere,
)

# ---------------------------------------------------------------------------
# Create primitives
# ---------------------------------------------------------------------------

sphere = implicit_sphere(center=(0, 0, 0), radius=0.45)
box = implicit_box(center=(0, 0, 0), half_extents=(0.3, 0.3, 0.3))

BOUNDS = (-0.8, 0.8) * 3
RES = 80

# ---------------------------------------------------------------------------
# Boolean operations
# ---------------------------------------------------------------------------

union = sphere | box
intersection = sphere & box
difference = sphere - box
smooth = sphere.smooth_union(box, k=0.3)

# ---------------------------------------------------------------------------
# Visualize in a 2x2 grid
# ---------------------------------------------------------------------------

plotter = pv.Plotter(shape=(2, 2), window_size=(1200, 1200))

ops = [
    (0, 0, "Union (A | B)", union, "steelblue"),
    (0, 1, "Intersection (A & B)", intersection, "coral"),
    (1, 0, "Difference (A - B)", difference, "mediumseagreen"),
    (1, 1, "Smooth Union (k=0.3)", smooth, "orchid"),
]

for row, col, title, shape, color in ops:
    mesh = shape.generate_vtk(bounds=BOUNDS, resolution=RES)
    plotter.subplot(row, col)
    plotter.add_text(title, font_size=12)
    plotter.add_mesh(mesh, color=color, show_edges=False)

plotter.link_views()
plotter.show()

# ---------------------------------------------------------------------------
# Batch smooth union — combining many primitives efficiently
# ---------------------------------------------------------------------------

import numpy as np

rng = np.random.default_rng(42)
spheres = [
    implicit_sphere(center=tuple(rng.uniform(-0.5, 0.5, 3)), radius=0.15)
    for _ in range(20)
]

combined = batch_smooth_union(spheres, k=0.1)
mesh = combined.generate_vtk(bounds=(-0.8, 0.8) * 3, resolution=100)

plotter2 = pv.Plotter()
plotter2.add_text("Batch smooth union (20 spheres, k=0.1)", font_size=12)
plotter2.add_mesh(mesh, color="coral", show_edges=False)
plotter2.show()
