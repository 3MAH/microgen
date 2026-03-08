"""Implicit (F-rep) lattice construction and comparison with strut lattices.

This example builds BCC, Octet-Truss, and Cubic lattices using the F-rep
implicit modeling framework (capsule primitives + smooth union) and compares
them side-by-side with the original strut-based lattice meshes from microgen.

A second comparison shows the effect of varying the smoothness parameter *k*
on a BCC lattice — from sharp joints (k ≈ 0) to heavily rounded blends.
"""

from itertools import product

import numpy as np
import pyvista as pv

from microgen import BodyCenteredCubic, Cubic, OctetTruss
from microgen.shape.implicit_shape import (
    ImplicitShape,
    batch_smooth_union,
    implicit_capsule,
    implicit_sphere,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

CELL = 1.0  # unit cell size
HALF = CELL / 2.0
STRUT_RADIUS = 0.08
RESOLUTION = 80
BOUNDS = (-HALF * 1.05, HALF * 1.05) * 3  # slightly larger than unit cube


def _cube_corners():
    return np.array(list(product([-HALF, HALF], repeat=3)))


def _face_centers():
    centers = []
    for axis in range(3):
        for sign in [-1, 1]:
            pt = [0.0, 0.0, 0.0]
            pt[axis] = sign * HALF
            centers.append(pt)
    return np.array(centers)


def _build_implicit_lattice(
    vertices: np.ndarray,
    pairs: list[tuple[int, int]],
    radius: float,
    k: float = 0.0,
    add_joints: bool = True,
) -> ImplicitShape:
    """Build an implicit lattice from vertices and connectivity.

    Each strut is an implicit capsule; joints are implicit spheres.
    All primitives are combined with ``batch_smooth_union`` to avoid
    recursion-depth issues with large primitive counts.
    """
    primitives: list[ImplicitShape] = []

    for i, j in pairs:
        primitives.append(
            implicit_capsule(
                start=tuple(vertices[i]),
                end=tuple(vertices[j]),
                radius=radius,
            )
        )

    if add_joints:
        for v in vertices:
            primitives.append(implicit_sphere(center=tuple(v), radius=radius))

    return batch_smooth_union(primitives, k=k)


# ---------------------------------------------------------------------------
# Lattice definitions (same topology as microgen strut lattices)
# ---------------------------------------------------------------------------


def cubic_lattice(
    radius: float = STRUT_RADIUS, k: float = 0.0,
) -> ImplicitShape:
    """Cubic lattice: 8 corner nodes, 12 edge struts."""
    verts = _cube_corners()
    pairs = []
    for i in range(len(verts)):
        for j in range(i + 1, len(verts)):
            if np.linalg.norm(verts[i] - verts[j]) < CELL + 1e-6:
                pairs.append((i, j))
    return _build_implicit_lattice(verts, pairs, radius, k)


def bcc_lattice(
    radius: float = STRUT_RADIUS, k: float = 0.0,
) -> ImplicitShape:
    """Body-centered cubic: center node connected to 8 corners."""
    corners = _cube_corners()
    verts = np.vstack(([[0.0, 0.0, 0.0]], corners))
    pairs = [(0, i) for i in range(1, 9)]
    return _build_implicit_lattice(verts, pairs, radius, k)


def octet_truss_lattice(
    radius: float = STRUT_RADIUS, k: float = 0.0,
) -> ImplicitShape:
    """Octet-truss: 8 corners + 6 face centers, connected at distance √2/2."""
    corners = _cube_corners()
    faces = _face_centers()
    verts = np.vstack((corners, faces))
    connection_dist = CELL / np.sqrt(2) + 1e-5
    pairs = []
    for i in range(len(verts)):
        for j in range(i + 1, len(verts)):
            if np.linalg.norm(verts[i] - verts[j]) < connection_dist:
                pairs.append((i, j))
    return _build_implicit_lattice(verts, pairs, radius, k)


# ---------------------------------------------------------------------------
# 1. Comparison: original strut lattice vs. implicit lattice (k=0 and k>0)
# ---------------------------------------------------------------------------

print("Building comparison: strut lattices vs implicit lattices ...")

lattice_defs = [
    ("Cubic", Cubic, cubic_lattice),
    ("BCC", BodyCenteredCubic, bcc_lattice),
    ("Octet-Truss", OctetTruss, octet_truss_lattice),
]

plotter = pv.Plotter(shape=(len(lattice_defs), 3), window_size=(1800, 1800))

for row, (name, StrutClass, implicit_fn) in enumerate(lattice_defs):
    # Column 0: original strut lattice
    strut = StrutClass(strut_radius=STRUT_RADIUS)
    strut_mesh = strut.generate_vtk()

    plotter.subplot(row, 0)
    plotter.add_text(f"{name}\n(strut)", font_size=10)
    plotter.add_mesh(strut_mesh, color="steelblue", show_edges=False)

    # Column 1: implicit lattice, hard union (k=0)
    imp_hard = implicit_fn(radius=STRUT_RADIUS, k=0.0)
    mesh_hard = imp_hard.generate_vtk(bounds=BOUNDS, resolution=RESOLUTION)

    plotter.subplot(row, 1)
    plotter.add_text(f"{name}\nimplicit k=0", font_size=10)
    plotter.add_mesh(mesh_hard, color="coral", show_edges=False)

    # Column 2: implicit lattice, smooth union (k=0.15)
    imp_smooth = implicit_fn(radius=STRUT_RADIUS, k=0.15)
    mesh_smooth = imp_smooth.generate_vtk(bounds=BOUNDS, resolution=RESOLUTION)

    plotter.subplot(row, 2)
    plotter.add_text(f"{name}\nimplicit k=0.15", font_size=10)
    plotter.add_mesh(mesh_smooth, color="mediumseagreen", show_edges=False)

plotter.link_views()
plotter.show()


# ---------------------------------------------------------------------------
# 2. Smoothness sweep on BCC lattice
# ---------------------------------------------------------------------------

print("Building smoothness sweep on BCC lattice ...")

k_values = [0.0, 0.05, 0.10, 0.20, 0.35]
n_k = len(k_values)

plotter2 = pv.Plotter(shape=(1, n_k), window_size=(400 * n_k, 500))

for col, k in enumerate(k_values):
    bcc = bcc_lattice(radius=STRUT_RADIUS, k=k)
    mesh = bcc.generate_vtk(bounds=BOUNDS, resolution=RESOLUTION)

    plotter2.subplot(0, col)
    plotter2.add_text(f"k = {k:.2f}", font_size=12)
    plotter2.add_mesh(mesh, color="coral", show_edges=False)

plotter2.link_views()
plotter2.show()
