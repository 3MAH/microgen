"""
Periodic RVE of short fibres
============================

A short fibre composite as a periodic RVE: randomly oriented fibres placed
without overlap in a unit cell, with their periodic images. Each fibre is an
implicit :class:`microgen.Capsule` (a cylinder with hemispherical ends) or a
:class:`microgen.Cylinder` (flat, sharp ends), and the reinforcement is the
union of their fields.

The same approach mixes inclusion types: spheres are added at the end, and
the hybrid RVE is again a single implicit shape.
"""

# sphinx_gallery_thumbnail_number = 3
import itertools

import numpy as np
import pyvista as pv
from scipy.spatial.transform import Rotation

from microgen import Capsule, Cylinder, Phase, Sphere, batch_smooth_union

# %%
# Random sequential addition of fibres of length ``length`` and radius
# ``radius``. A fibre is accepted if the points sampled along its axis stay
# at more than ``2 * radius + gap`` from the points of the accepted fibres,
# with the minimum image convention across the faces.
length, radius, gap = 0.4, 0.04, 0.01
n_fibres = 30

rng = np.random.default_rng(2)
t = np.linspace(-0.5, 0.5, 15)[:, None]  # sampling along the axis
fibres = []  # (centre, Euler angles, axis points)
accepted_points = np.empty((0, 3))
while len(fibres) < n_fibres:
    center = rng.uniform(0.0, 1.0, 3)
    rotation = Rotation.random(random_state=rng)
    axis = rotation.apply((1.0, 0.0, 0.0))
    points = center + length * t * axis
    if len(accepted_points):
        delta = points[:, None, :] - accepted_points[None, :, :]
        delta -= np.round(delta)
        if np.linalg.norm(delta, axis=2).min() <= 2 * radius + gap:
            continue
    fibres.append((center, rotation.as_euler("ZXZ", degrees=True)))
    accepted_points = np.vstack([accepted_points, points])

# %%
# One fibre shape per fibre and per periodic image that reaches the cell.
reach = 0.5 * length + radius
images = [
    (tuple(center + shift), tuple(angles))
    for center, angles in fibres
    for shift in itertools.product((-1, 0, 1), repeat=3)
    if np.all((center + shift > -reach) & (center + shift < 1 + reach))
]
print(f"{n_fibres} fibres, {len(images)} with their periodic images")

# %%
# Rounded fibres
# --------------
# Capsules: the axis of length ``length`` is capped by two hemispheres.
capsules = [
    Capsule(center=center, orientation=angles, height=length, radius=radius)
    for center, angles in images
]

# The phase is sampled on a grid whose nodes lie on the faces of the cell,
# and its surface is taken from the volume cells where the field is
# negative: fibres cut by a face get flat caps with sharp edges.
cell_bounds = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
fibre_phase = Phase.from_shape(batch_smooth_union(capsules), bounds=cell_bounds, resolution=121)


def surface_in_cell(phase):
    """Surface of the cells where the field is negative, cut flat by the cell faces."""
    volume = phase.volume_mesh()
    return volume.extract_surface(algorithm="dataset_surface").triangulate().clean()


single = np.pi * radius**2 * length + 4 / 3 * np.pi * radius**3
measured = fibre_phase.volume_mesh().volume
print(f"Fibre volume fraction: {measured:.4f} (exact {n_fibres * single:.4f})")
fibre_surface = surface_in_cell(fibre_phase)
fibre_surface.save("fibres_rve.stl")

box = pv.Box(bounds=(0, 1, 0, 1, 0, 1))
pl = pv.Plotter()
pl.add_mesh(fibre_surface, color="#ef7d2d", smooth_shading=True, split_sharp_edges=True)
pl.add_mesh(box, style="wireframe", color="black")
pl.show()

# %%
# Sharp fibres
# ------------
# Cylinders with flat ends, on the same axes. A cylinder of height ``length``
# lies inside the capsule of the same fibre, so the fibres still do not
# overlap. ``split_sharp_edges`` keeps the edges of the ends sharp in the
# rendering.
cylinders = [
    Cylinder(center=center, orientation=angles, height=length, radius=radius)
    for center, angles in images
]
sharp_phase = Phase.from_shape(batch_smooth_union(cylinders), bounds=cell_bounds, resolution=121)
measured = sharp_phase.volume_mesh().volume
exact = n_fibres * np.pi * radius**2 * length
print(f"Sharp fibre volume fraction: {measured:.4f} (exact {exact:.4f})")

pl = pv.Plotter()
pl.add_mesh(
    surface_in_cell(sharp_phase), color="#ef7d2d", smooth_shading=True, split_sharp_edges=True
)
pl.add_mesh(box, style="wireframe", color="black")
pl.show()

# %%
# Hybrid RVE
# ----------
# Add spheres in the space left by the fibres. The test uses the
# fibre field itself: a sphere is accepted if the field of the fibres at its
# centre (the distance to the nearest fibre surface) exceeds its radius.
fibre_field = batch_smooth_union(capsules)
sphere_radius = 0.07
sphere_centers = []
while len(sphere_centers) < 12:
    candidate = rng.uniform(0.0, 1.0, 3)
    clear_of_fibres = fibre_field.evaluate(*candidate[:, None]) > sphere_radius + gap
    delta = np.asarray(sphere_centers).reshape(-1, 3) - candidate
    delta -= np.round(delta)
    clear_of_spheres = np.all(np.linalg.norm(delta, axis=1) > 2 * sphere_radius + gap)
    if clear_of_fibres[0] and clear_of_spheres:
        sphere_centers.append(candidate)

spheres = [
    Sphere(center=tuple(center + shift), radius=sphere_radius)
    for center in sphere_centers
    for shift in itertools.product((-1, 0, 1), repeat=3)
    if np.all((center + shift > -sphere_radius) & (center + shift < 1 + sphere_radius))
]
particle_phase = Phase.from_shape(batch_smooth_union(spheres), bounds=cell_bounds, resolution=121)

# %%
# The hybrid RVE: fibres in orange, spheres in blue.
pl = pv.Plotter()
pl.add_mesh(fibre_surface, color="#ef7d2d", smooth_shading=True, split_sharp_edges=True)
pl.add_mesh(
    surface_in_cell(particle_phase), color="#6a8caf", smooth_shading=True, split_sharp_edges=True
)
pl.add_mesh(box, style="wireframe", color="black")
pl.show()
