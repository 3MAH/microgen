"""Generate a trabecular bone-like structure from Voronoi edges.

Seed points are scattered in a cube, the 3D Voronoi tessellation provides
the connectivity (edges with natural degree-4 nodes), and each edge becomes
an implicit capsule.  All primitives are combined with smooth union and
intersected with a box for a clean cross-section.

Reference: Voronoi tessellation is the standard approach in literature for
synthetic trabecular bone scaffolds — interior Voronoi vertices have exactly
4 edges, matching real trabecular node valence.
"""

import numpy as np
import pyvista as pv
from scipy.spatial import Voronoi

from microgen.shape.implicit_shape import (
    batch_smooth_union,
    implicit_box,
    implicit_capsule,
    implicit_sphere,
)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

SEED = 42
N_POINTS = 150  # Voronoi seed count

# Domain
DOMAIN_SIZE = 2.0
HALF = DOMAIN_SIZE / 2.0

# Clip box — the clean cut
BOX_HALF = 0.75

# Strut geometry (trabecular dimensions)
RADIUS_MEAN = 0.035  # ~Tb.Th / 2
RADIUS_STD = 0.006

# Smooth blending
K = 0.03

# Meshing
RESOLUTION = 200
BOUNDS = (-BOX_HALF * 1.02, BOX_HALF * 1.02) * 3

# Lloyd relaxation iterations (regularizes cell sizes while keeping randomness)
LLOYD_ITERATIONS = 2

# ---------------------------------------------------------------------------
# Generate seed points + Lloyd relaxation
# ---------------------------------------------------------------------------

rng = np.random.default_rng(SEED)

# Scatter seeds in a region larger than the clip box
seeds = rng.uniform(-HALF, HALF, size=(N_POINTS, 3))


def lloyd_relaxation(pts, bounds_min, bounds_max, iterations=1):
    """Move each seed toward the centroid of its Voronoi cell."""
    for _ in range(iterations):
        vor = Voronoi(pts)
        new_pts = np.copy(pts)
        for i, region_idx in enumerate(vor.point_region):
            region = vor.regions[region_idx]
            if -1 in region or len(region) == 0:
                continue  # skip open/infinite regions
            verts = vor.vertices[region]
            centroid = verts.mean(axis=0)
            centroid = np.clip(centroid, bounds_min, bounds_max)
            new_pts[i] = centroid
        pts = new_pts
    return pts


print(f"Lloyd relaxation ({LLOYD_ITERATIONS} iterations) on {N_POINTS} seeds ...")
seeds = lloyd_relaxation(seeds, -HALF, HALF, iterations=LLOYD_ITERATIONS)

# ---------------------------------------------------------------------------
# Compute Voronoi and extract edges inside the clip box
# ---------------------------------------------------------------------------

print("Computing Voronoi tessellation ...")
vor = Voronoi(seeds)

margin = BOX_HALF + 0.15
edge_set = set()
for ridge in vor.ridge_vertices:
    if -1 in ridge:
        continue
    n = len(ridge)
    for i in range(n):
        v1, v2 = ridge[i], ridge[(i + 1) % n]
        edge_set.add((min(v1, v2), max(v1, v2)))

# Filter: keep edges near the clip box
edges = []
for v1, v2 in edge_set:
    p1 = vor.vertices[v1]
    p2 = vor.vertices[v2]
    if np.all(np.abs(p1) < margin) or np.all(np.abs(p2) < margin):
        edges.append((v1, v2))

print(f"  {len(edges)} finite edges near clip box")

# Collect unique vertex indices used
used_verts = set()
for v1, v2 in edges:
    used_verts.add(v1)
    used_verts.add(v2)

print(f"  {len(used_verts)} Voronoi vertices (nodes)")

# ---------------------------------------------------------------------------
# Build implicit primitives
# ---------------------------------------------------------------------------

print(f"Building {len(edges)} capsules + {len(used_verts)} joint spheres ...")

radii = np.clip(rng.normal(RADIUS_MEAN, RADIUS_STD, size=len(edges)), 0.015, None)

primitives = []

for idx, (v1, v2) in enumerate(edges):
    primitives.append(
        implicit_capsule(
            start=tuple(vor.vertices[v1]),
            end=tuple(vor.vertices[v2]),
            radius=float(radii[idx]),
        )
    )

for vi in used_verts:
    primitives.append(
        implicit_sphere(center=tuple(vor.vertices[vi]), radius=RADIUS_MEAN * 1.1)
    )

# ---------------------------------------------------------------------------
# Combine with smooth union + box intersection
# ---------------------------------------------------------------------------

print(f"Combining {len(primitives)} primitives with batch smooth union (k={K}) ...")

result = batch_smooth_union(primitives, k=K)

print(f"Intersecting with box (half-extent={BOX_HALF}) ...")
box = implicit_box(center=(0, 0, 0), half_extents=(BOX_HALF, BOX_HALF, BOX_HALF))
result = result & box

# ---------------------------------------------------------------------------
# Mesh and visualize
# ---------------------------------------------------------------------------

print(f"Meshing (resolution={RESOLUTION}) ...")
mesh = result.generate_vtk(bounds=BOUNDS, resolution=RESOLUTION)

plotter = pv.Plotter()
plotter.add_text("Trabecular Bone (Voronoi)", font_size=14)
plotter.add_mesh(mesh, color="ivory", show_edges=False)
plotter.show()
