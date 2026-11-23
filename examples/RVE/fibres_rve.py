"""
Periodic RVE of short fibres
============================

A short fibre composite as a representative volume element (RVE): randomly
oriented fibres placed without overlap in a unit cell. Each fibre is an
implicit :class:`microgen.Capsule` (a cylinder with hemispherical ends) or a
:class:`microgen.Cylinder` (flat, sharp ends), and the reinforcement is the
union of their fields.

The switch ``periodic`` below selects the kind of RVE:

- ``periodic = True`` (default): the microstructure is periodic. A fibre that
  crosses a face of the cell reappears on the opposite face, so that copies
  of the RVE tile space without a seam. This is the RVE needed for periodic
  boundary conditions, as in the homogenization tools of fedoo and simcoon.
- ``periodic = False``: the fibres are simply cut by the faces of the cell
  and do not reappear on the other side. Opposite faces do not match, and
  the RVE is meant for other boundary conditions (imposed displacements or
  tractions on each face).

Periodicity takes two ingredients, both switched by ``periodic``: distances
between fibres measured across the faces (minimum image convention) when the
fibres are placed, and periodic images of the fibres when they are built.

The same approach mixes inclusion types: spheres are added at the end, and
the hybrid RVE is again a single implicit shape.
"""

# sphinx_gallery_thumbnail_number = 4
import itertools

import numpy as np
import pyvista as pv
from scipy.spatial.transform import Rotation

from microgen import Capsule, Cylinder, Phase, Sphere, batch_smooth_union

# %%
# Parameters
# ----------
# ``periodic`` selects a periodic RVE (``True``) or fibres cut by the cell
# (``False``); every step below that depends on it says so.
periodic = True
length, radius, gap = 0.4, 0.04, 0.01
n_fibres = 30

# %%
# Placing the fibres
# ------------------
# Random sequential addition: a fibre is drawn at random and accepted if the
# points sampled along its axis stay at more than ``2 * radius + gap`` from
# the points of the accepted fibres.
#
# *Periodic*: the distance is measured with the minimum image convention,
# ``delta -= round(delta)``, i.e. to the nearest copy of the other fibre in
# the neighbouring cells. Two fibres near opposite faces are then also kept
# apart, since they will meet once the fibres wrap around.
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
        if periodic:
            delta -= np.round(delta)
        if np.linalg.norm(delta, axis=2).min() <= 2 * radius + gap:
            continue
    fibres.append((center, rotation.as_euler("ZXZ", degrees=True)))
    accepted_points = np.vstack([accepted_points, points])

# %%
# Periodic images
# ---------------
# *Periodic*: each fibre is copied with the shifts of the 26 neighbouring
# cells, ``(i, j, k)`` with ``i, j, k`` in ``{-1, 0, 1}``, and the copies
# that reach the cell are kept. The part of a fibre that leaves the cell
# through the face x = 1 is thus brought back through the face x = 0 by its
# copy shifted by ``(-1, 0, 0)``.
#
# *Not periodic*: only the fibres themselves are kept, ``shift = (0, 0, 0)``.
shifts = list(itertools.product((-1, 0, 1), repeat=3)) if periodic else [(0, 0, 0)]
reach = 0.5 * length + radius
images = [
    (tuple(center + shift), tuple(angles))
    for center, angles in fibres
    for shift in shifts
    if np.all((center + shift > -reach) & (center + shift < 1 + reach))
]
print(f"{n_fibres} fibres, {len(images)} shapes with their periodic images")

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
fibre_phase = Phase.from_shape(
    batch_smooth_union(capsules), bounds=cell_bounds, resolution=121
)


def surface_in_cell(phase):
    """Surface of the cells where the field is negative, cut flat by the cell faces."""
    volume = phase.volume_mesh()
    return volume.extract_surface(algorithm="dataset_surface").triangulate().clean()


# *Periodic*: the parts of a fibre cut by a face come back through the
# opposite face, so the cell holds exactly ``n_fibres`` whole fibres and the
# volume fraction is known in closed form. *Not periodic*: the parts outside
# the cell are lost and the fraction is lower.
single = np.pi * radius**2 * length + 4 / 3 * np.pi * radius**3
measured = fibre_phase.volume_mesh().volume
if periodic:
    print(f"Fibre volume fraction: {measured:.4f} (exact {n_fibres * single:.4f})")
else:
    print(
        f"Fibre volume fraction: {measured:.4f} (whole fibres: {n_fibres * single:.4f})"
    )
fibre_surface = surface_in_cell(fibre_phase)
fibre_surface.save("fibres_rve.stl")

box = pv.Box(bounds=(0, 1, 0, 1, 0, 1))
pl = pv.Plotter()
pl.add_mesh(fibre_surface, color="#ef7d2d", smooth_shading=True, split_sharp_edges=True)
pl.add_mesh(box, style="wireframe", color="black")
pl.show()

# %%
# Checking the periodicity
# ------------------------
# The implicit field gives a direct check: a point of a face is in a fibre
# exactly when its image on the opposite face is. The test holds on the
# three pairs of faces of a periodic RVE, and fails when ``periodic = False``.
fibre_field = batch_smooth_union(capsules)
u, v = (a.ravel() for a in np.meshgrid(np.linspace(0, 1, 80), np.linspace(0, 1, 80)))
zero, one = np.zeros_like(u), np.ones_like(u)
for axis, face_0, face_1 in (
    ("x", (zero, u, v), (one, u, v)),
    ("y", (u, zero, v), (u, one, v)),
    ("z", (u, v, zero), (u, v, one)),
):
    same = np.array_equal(
        fibre_field.evaluate(*face_0) < 0, fibre_field.evaluate(*face_1) < 0
    )
    print(f"Faces {axis} = 0 and {axis} = 1 match: {same}")

# %%
# Sharp fibres
# ------------
# Cylinders with flat ends, on the same axes. A cylinder of height ``length``
# lies inside the capsule of the same fibre, so the fibres still do not
# overlap.
#
# Marching cubes round off every sharp edge (the flat ends, and the cuts by
# the cell faces). Since the fibres do not overlap, each one can instead be
# meshed exactly with :meth:`microgen.Cylinder.generate_surface_mesh` and
# clipped by the six faces of the cell, which keeps every edge sharp.
cylinders = [
    Cylinder(center=center, orientation=angles, height=length, radius=radius)
    for center, angles in images
]
cell_faces = [
    ((-1, 0, 0), (1, 0, 0)),
    ((1, 0, 0), (0, 0, 0)),
    ((0, -1, 0), (0, 1, 0)),
    ((0, 1, 0), (0, 0, 0)),
    ((0, 0, -1), (0, 0, 1)),
    ((0, 0, 1), (0, 0, 0)),
]  # (normal pointing into the cell, point on the face)
pieces = []
for cylinder in cylinders:
    mesh = cylinder.generate_surface_mesh(resolution=48).triangulate().clean()
    for normal, origin in cell_faces:
        # Signed distance of the points to the face, positive inside the cell.
        inside = (mesh.points - origin) @ np.asarray(normal, dtype=float)
        if inside.max() <= 0.0:  # fibre image entirely out of the cell
            mesh = None
            break
        if inside.min() < 0.0:  # fibre crossing the face
            mesh = mesh.clip_closed_surface(normal=normal, origin=origin).clean()
    if mesh is not None:
        pieces.append(mesh)
sharp_fibres = pv.merge(pieces)
sharp_fibres.save("fibres_rve_sharp.stl")

exact = n_fibres * np.pi * radius**2 * length
label = "exact" if periodic else "whole fibres"
print(f"Sharp fibre volume fraction: {sharp_fibres.volume:.4f} ({label} {exact:.4f})")

pl = pv.Plotter()
pl.add_mesh(sharp_fibres, color="#ef7d2d", smooth_shading=True, split_sharp_edges=True)
pl.add_mesh(box, style="wireframe", color="black")
pl.show()

# %%
# Tiling the RVE
# --------------
# Four copies of the RVE side by side. *Periodic*: the fibres cut by a face
# continue in the neighbouring cell and the tiling shows no seam. *Not
# periodic*: the fibres stop at the faces and the seams are visible.
pl = pv.Plotter()
for i, j in itertools.product(range(2), repeat=2):
    shift = (float(i), float(j), 0.0)
    pl.add_mesh(
        sharp_fibres.translate(shift, inplace=False),
        color="#ef7d2d",
        smooth_shading=True,
        split_sharp_edges=True,
    )
    pl.add_mesh(box.translate(shift, inplace=False), style="wireframe", color="black")
pl.show()

# %%
# Hybrid RVE
# ----------
# Add spheres in the space left by the fibres. The test uses the
# fibre field itself: a sphere is accepted if the field of the fibres at its
# centre (the distance to the nearest fibre surface) exceeds its radius. The
# spheres follow the same ``periodic`` switch: minimum image distance between
# spheres, and periodic images.
sphere_radius = 0.07
sphere_centers = []
while len(sphere_centers) < 12:
    candidate = rng.uniform(0.0, 1.0, 3)
    clear_of_fibres = fibre_field.evaluate(*candidate[:, None]) > sphere_radius + gap
    delta = np.asarray(sphere_centers).reshape(-1, 3) - candidate
    if periodic:
        delta -= np.round(delta)
    clear_of_spheres = np.all(np.linalg.norm(delta, axis=1) > 2 * sphere_radius + gap)
    if clear_of_fibres[0] and clear_of_spheres:
        sphere_centers.append(candidate)

spheres = [
    Sphere(center=tuple(center + shift), radius=sphere_radius)
    for center in sphere_centers
    for shift in shifts
    if np.all((center + shift > -sphere_radius) & (center + shift < 1 + sphere_radius))
]
particle_phase = Phase.from_shape(
    batch_smooth_union(spheres), bounds=cell_bounds, resolution=121
)

# %%
# The hybrid RVE: fibres in orange, spheres in blue.
pl = pv.Plotter()
pl.add_mesh(fibre_surface, color="#ef7d2d", smooth_shading=True, split_sharp_edges=True)
pl.add_mesh(
    surface_in_cell(particle_phase),
    color="#6a8caf",
    smooth_shading=True,
    split_sharp_edges=True,
)
pl.add_mesh(box, style="wireframe", color="black")
pl.show()
