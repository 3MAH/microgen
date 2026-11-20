"""
Spherical, cylindrical and swept frames
=======================================

Three geometries built with the parametric-grid TPMS classes, each clipped by
a :math:`y` plane so the interior is visible:

* :class:`~microgen.shape.tpms.SphericalTpms`: a gyroid wrapping a sphere
  (``sphere_g4_tpms.vtk``);
* :class:`~microgen.shape.tpms.CylindricalTpms`: a gyroid wrapping a cylinder
  (``cylinder_g5_tpms.vtk``);
* :class:`~microgen.shape.tpms.Sweep`: a gyroid following a helical curve
  (``sweep_g7_helix.vtk``).

Parameters:

- ``CELL_SIZE``        — uniform cubic cell edge (scalar → no stretched cells)
- ``OFFSET``           — sheet thickness
- ``SPHERE_RADIUS``    — radius of the spherical demo
- ``CYLINDER_RADIUS``, ``CYLINDER_HEIGHT`` — for the cylindrical demo
"""

from __future__ import annotations

import numpy as np
import pyvista as pv

from microgen.shape.surface_functions import gyroid
from microgen.shape.tpms import CylindricalTpms, SphericalTpms, Sweep


# -----------------------------------------------------------------------------
# Knobs
# -----------------------------------------------------------------------------

CELL_SIZE = 1.0
OFFSET = 0.5
SPHERE_RADIUS = 3.0
CYLINDER_RADIUS = 1.5
CYLINDER_HEIGHT = 6.0


def clip_y(mesh: pv.DataSet, frac: float = -0.2) -> pv.DataSet:
    """Cut a mesh by a y-plane at ``frac × bbox_y`` so internals are visible."""
    bb = np.array(mesh.bounds)
    y_origin = bb[2] + frac * (bb[3] - bb[2])
    return mesh.clip("y", origin=(0.0, y_origin, 0.0), invert=False)


def report(label: str, mesh: pv.DataSet) -> None:
    print(f"{label:32s}  vol={abs(mesh.volume):7.3f}  n_cells={mesh.n_cells}")


# -----------------------------------------------------------------------------
# 1. SphericalTpms — gyroid wrapping a sphere (auto-fill θ + φ via repeat=0)
# -----------------------------------------------------------------------------

sph = SphericalTpms(
    radius=SPHERE_RADIUS,
    surface_function=gyroid,
    offset=OFFSET,
    cell_size=CELL_SIZE,
    repeat_cell=(2, 0, 0),
)
m_sph = sph.generate_surface_mesh(type_part="sheet")
clip_y(m_sph).save("sphere_g4_tpms.vtk")
report("1. SphericalTpms (R=3)", m_sph)

# -----------------------------------------------------------------------------
# 2. CylindricalTpms — gyroid wrapping a cylinder (auto-fill θ via repeat=0)
# -----------------------------------------------------------------------------

cyl = CylindricalTpms(
    radius=CYLINDER_RADIUS,
    surface_function=gyroid,
    offset=OFFSET,
    cell_size=CELL_SIZE,
    repeat_cell=(2, 0, int(CYLINDER_HEIGHT / CELL_SIZE)),
)
m_cyl = cyl.generate_surface_mesh(type_part="sheet")
clip_y(m_cyl).save("cylinder_g5_tpms.vtk")
report("2. CylindricalTpms (R=1.5)", m_cyl)

# -----------------------------------------------------------------------------
# 3. Sweep — gyroid along a helical curve (1.5 turns, height 6)
# -----------------------------------------------------------------------------


def helix(t: float) -> np.ndarray:
    """Helix: 1.5 turns, radius 2, height 6."""
    theta = 2.0 * np.pi * 1.5 * t
    return np.array([2.0 * np.cos(theta), 2.0 * np.sin(theta), 6.0 * (t - 0.5)])


sw = Sweep(
    curve_points=helix,
    surface_function=gyroid,
    radial_max=0.6,
    offset=OFFSET,
    cell_size=CELL_SIZE,
    repeat_cell=(8, 1, 6),
    n_curve_samples=200,
)
m_sw = sw.generate_surface_mesh(type_part="sheet")
clip_y(m_sw).save("sweep_g7_helix.vtk")
report("3. Sweep along helix", m_sw)

# %%
# Display the three clipped geometries side by side.
pl = pv.Plotter(shape=(1, 3), window_size=(1500, 500))
for i, mesh in enumerate((m_sph, m_cyl, m_sw)):
    pl.subplot(0, i)
    pl.add_mesh(clip_y(mesh), smooth_shading=True, color="#ef7d2d")
    pl.view_xz()
pl.link_views()
pl.show()
