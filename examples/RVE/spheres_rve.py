"""
Periodic RVE of spheres
=======================

A particle-reinforced composite as a representative volume element (RVE):
spheres placed at random without overlap in a unit cell. Every sphere is an
implicit :class:`microgen.Sphere`; their union is one minimum over the
fields (:func:`microgen.batch_smooth_union`) and the matrix is its
complement. No CAD kernel is involved.

The switch ``periodic`` below selects the kind of RVE:

- ``periodic = True`` (default): the microstructure is periodic. A sphere
  that crosses a face of the cell reappears on the opposite face, so that
  copies of the RVE tile space without a seam. This is the RVE needed for
  periodic boundary conditions, as in the homogenization tools of fedoo and
  simcoon.
- ``periodic = False``: the spheres are simply cut by the faces of the cell
  and do not reappear on the other side. Opposite faces do not match, and
  the RVE is meant for other boundary conditions (imposed displacements or
  tractions on each face).

Periodicity takes two ingredients, both switched by ``periodic``: distances
between spheres measured across the faces (minimum image convention) when
the spheres are placed, and periodic images of the spheres when they are
built.
"""

# sphinx_gallery_thumbnail_number = 1
import itertools

import numpy as np
import pyvista as pv

from microgen import Phase, Sphere, batch_smooth_union

# %%
# Parameters
# ----------
# ``periodic`` selects a periodic RVE (``True``) or spheres cut by the cell
# (``False``); every step below that depends on it says so.
periodic = True
radius = 0.12
gap = 0.02
n_spheres = 25

# %%
# Placing the spheres
# -------------------
# Random sequential addition: draw centres one by one and keep a centre if it
# is at more than ``2 * radius + gap`` from the accepted ones.
#
# *Periodic*: the distance is measured with the minimum image convention,
# ``delta -= round(delta)``, i.e. to the nearest copy of the other sphere in
# the neighbouring cells. Two spheres near opposite faces are then also kept
# apart, since they will meet once the spheres wrap around.
rng = np.random.default_rng(0)
centers = []
while len(centers) < n_spheres:
    candidate = rng.uniform(0.0, 1.0, 3)
    delta = np.asarray(centers).reshape(-1, 3) - candidate
    if periodic:
        delta -= np.round(delta)
    if np.all(np.linalg.norm(delta, axis=1) > 2 * radius + gap):
        centers.append(candidate)

# %%
# Periodic images
# ---------------
# *Periodic*: each sphere is copied with the shifts of the 26 neighbouring
# cells, ``(i, j, k)`` with ``i, j, k`` in ``{-1, 0, 1}``, and the copies
# that reach the cell are kept. The part of a sphere that leaves the cell
# through the face x = 1 is thus brought back through the face x = 0 by its
# copy shifted by ``(-1, 0, 0)``.
#
# *Not periodic*: only the spheres themselves are kept, ``shift = (0, 0, 0)``.
shifts = list(itertools.product((-1, 0, 1), repeat=3)) if periodic else [(0, 0, 0)]
spheres = [
    Sphere(center=tuple(center + np.array(shift)), radius=radius)
    for center in centers
    for shift in shifts
    if np.all((center + shift > -radius) & (center + shift < 1 + radius))
]
print(f"{n_spheres} spheres, {len(spheres)} shapes with their periodic images")

# %%
# Phases
# ------
# The inclusions are the union of the spheres; the matrix is its complement
# ``~inclusions`` (the boolean operators of :class:`microgen.Shape` work on
# the fields). Both phases are sampled on a grid whose nodes lie on the faces
# of the cell.
inclusions = batch_smooth_union(spheres)
cell_bounds = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
particles = Phase.from_shape(inclusions, bounds=cell_bounds, resolution=101)
matrix_phase = Phase.from_shape(~inclusions, bounds=cell_bounds, resolution=101)

# *Periodic*: the parts of a sphere cut by a face come back through the
# opposite face, so the cell holds exactly ``n_spheres`` whole spheres and
# the volume fraction is known in closed form. *Not periodic*: the parts
# outside the cell are lost and the fraction is lower.
whole = n_spheres * 4 / 3 * np.pi * radius**3
measured = particles.volume_mesh().volume
label = "exact" if periodic else "whole spheres"
print(f"Volume fraction: {measured:.4f} ({label} {whole:.4f})")

# %%
# Checking the periodicity
# ------------------------
# The implicit field gives a direct check: a point of a face is in a particle
# exactly when its image on the opposite face is. The test holds on the
# three pairs of faces of a periodic RVE, and fails when ``periodic = False``.
u, v = (a.ravel() for a in np.meshgrid(np.linspace(0, 1, 80), np.linspace(0, 1, 80)))
zero, one = np.zeros_like(u), np.ones_like(u)
for axis, face_0, face_1 in (
    ("x", (zero, u, v), (one, u, v)),
    ("y", (u, zero, v), (u, one, v)),
    ("z", (u, v, zero), (u, v, one)),
):
    same = np.array_equal(
        inclusions.evaluate(*face_0) < 0, inclusions.evaluate(*face_1) < 0
    )
    print(f"Faces {axis} = 0 and {axis} = 1 match: {same}")


# %%
# Surfaces
# --------
# The surfaces come from the volume cells where the field is negative: a
# particle cut by a face gets an exactly flat cap with a sharp edge, where a
# marching-cubes surface would round the edge off.
def surface_in_cell(phase):
    """Surface of the cells where the field is negative, cut flat by the cell faces."""
    volume = phase.volume_mesh()
    return volume.extract_surface(algorithm="dataset_surface").triangulate().clean()


particle_surface = surface_in_cell(particles)
particle_surface.save("spheres_rve.stl")

# %%
# The particles in the cell.
box = pv.Box(bounds=cell_bounds)
pl = pv.Plotter()
pl.add_mesh(
    particle_surface, color="#ef7d2d", smooth_shading=True, split_sharp_edges=True
)
pl.add_mesh(box, style="wireframe", color="black")
pl.show()

# %%
# The matrix, with the holes left by the particles on its faces.
surface_in_cell(matrix_phase).plot(
    color="#6a8caf", smooth_shading=True, split_sharp_edges=True
)

# %%
# Tiling the RVE
# --------------
# Four copies of the RVE side by side. *Periodic*: the spheres cut by a face
# continue in the neighbouring cell and the tiling shows no seam. *Not
# periodic*: the spheres stop at the faces and the seams are visible.
pl = pv.Plotter()
for i, j in itertools.product(range(2), repeat=2):
    shift = (float(i), float(j), 0.0)
    pl.add_mesh(
        particle_surface.translate(shift, inplace=False),
        color="#ef7d2d",
        smooth_shading=True,
        split_sharp_edges=True,
    )
    pl.add_mesh(box.translate(shift, inplace=False), style="wireframe", color="black")
pl.show()
