"""
Periodic RVE of spheres
=======================

A particle-reinforced composite as a periodic representative volume
element (RVE): spheres placed at random without overlap in a unit cell, with
their periodic images, so that a sphere cut by a face reappears on the
opposite face.

Every sphere is an implicit :class:`microgen.Sphere`; their union is one
minimum over the fields (:func:`microgen.batch_smooth_union`) and the matrix
is its complement. No CAD kernel is involved.
"""

# sphinx_gallery_thumbnail_number = 1
import itertools

import numpy as np
import pyvista as pv

from microgen import Phase, Sphere, batch_smooth_union

# %%
# Random sequential addition: draw centres one by one and keep a centre if it
# is far enough from the previous ones. Distances use the minimum image
# convention, so spheres do not overlap across the faces either.
radius = 0.12
gap = 0.02
n_spheres = 25

rng = np.random.default_rng(0)
centers = []
while len(centers) < n_spheres:
    candidate = rng.uniform(0.0, 1.0, 3)
    delta = np.asarray(centers).reshape(-1, 3) - candidate
    delta -= np.round(delta)  # minimum image in the unit cell
    if np.all(np.linalg.norm(delta, axis=1) > 2 * radius + gap):
        centers.append(candidate)

# %%
# Add the periodic images of the spheres that cross a face of the cell.
spheres = []
for center in centers:
    for shift in itertools.product((-1, 0, 1), repeat=3):
        image = center + np.array(shift)
        if np.all((image > -radius) & (image < 1 + radius)):
            spheres.append(Sphere(center=tuple(image), radius=radius))
print(f"{n_spheres} spheres, {len(spheres)} with their periodic images")

# %%
# The inclusions are the union of the spheres; the matrix is its complement
# ``~inclusions`` (the boolean operators of :class:`microgen.Shape` work on
# the fields). Both phases are sampled on a grid whose nodes lie on the faces
# of the cell.
inclusions = batch_smooth_union(spheres)
cell_bounds = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
particles = Phase.from_shape(inclusions, bounds=cell_bounds, resolution=101)
matrix_phase = Phase.from_shape(~inclusions, bounds=cell_bounds, resolution=101)


def surface_in_cell(phase):
    """Surface of the cells where the field is negative, cut flat by the cell faces."""
    volume = phase.volume_mesh()
    return volume.extract_surface(algorithm="dataset_surface").triangulate().clean()


target = n_spheres * 4 / 3 * np.pi * radius**3
measured = particles.volume_mesh().volume
print(f"Volume fraction: {measured:.4f} (exact {target:.4f})")

# %%
# A point of the x = 0 face is in a particle exactly when its image on the
# x = 1 face is.
y, z = np.meshgrid(np.linspace(0, 1, 50), np.linspace(0, 1, 50))
left = inclusions.evaluate(np.zeros(y.size), y.ravel(), z.ravel())
right = inclusions.evaluate(np.ones(y.size), y.ravel(), z.ravel())
print(f"Same phase on both faces: {np.array_equal(left < 0, right < 0)}")

# %%
# The surfaces come from the volume cells where the field is negative: a
# particle cut by a face gets an exactly flat cap with a sharp edge, where a
# marching-cubes surface would round the edge off.
particle_surface = surface_in_cell(particles)
particle_surface.save("spheres_rve.stl")

# %%
# The particles in the cell.
pl = pv.Plotter()
pl.add_mesh(
    particle_surface, color="#ef7d2d", smooth_shading=True, split_sharp_edges=True
)
pl.add_mesh(pv.Box(bounds=(0, 1, 0, 1, 0, 1)), style="wireframe", color="black")
pl.show()

# %%
# The matrix, with the holes left by the particles on its faces.
surface_in_cell(matrix_phase).plot(
    color="#6a8caf", smooth_shading=True, split_sharp_edges=True
)
