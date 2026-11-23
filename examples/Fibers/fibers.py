"""
Short fibre composite
=====================

Rebuild the fibres of a short fibre composite from a tomography analysis:
a CSV file gives the centre, diameter and orientation of about 15 000
fibres. Each fibre is a :class:`microgen.Cylinder`, whose orientation gives
the fibre axis; the fibres are drawn as PyVista tubes.

With the CAD backend (``pip install microgen[cad]``), the same cylinders can
be exported to STEP with ``Cylinder.generate_cad()`` and
:func:`microgen.cad.make_compound`, which takes several minutes for this
number of fibres.
"""

import numpy as np
import pyvista as pv

from microgen import Cylinder

# %%
# Read the fibres. Columns: centre (x, y, z), equivalent diameter and the
# two orientation angles, in µm and degrees.
data = np.loadtxt("fibers.csv", delimiter=",", skiprows=2)
centers = data[:, 0:3]
diameters = data[:, 3]
phi, theta = data[:, 4], data[:, 5]
print(f"{len(data)} fibres in a {np.ptp(centers, axis=0).round()} µm box")

# %%
# One :class:`~microgen.Cylinder` per fibre, of constant length. Its axis is
# the local x axis, rotated by the Euler angles (ZXZ) of the analysis.
length = 250.0  # µm
fibres = [
    Cylinder(
        center=tuple(c),
        orientation=(p, 90.0, t),
        height=length,
        radius=0.5 * d,
    )
    for c, d, p, t in zip(centers, diameters, phi, theta)
]
axes = np.array([fibre.orientation.apply((1.0, 0.0, 0.0)) for fibre in fibres])

# %%
# Draw each fibre as a tube between its two ends, coloured by the angle
# between the fibre and the z axis.
n = len(fibres)
ends = np.empty((2 * n, 3))
ends[0::2] = centers - 0.5 * length * axes
ends[1::2] = centers + 0.5 * length * axes
lines = np.column_stack([np.full(n, 2), np.arange(0, 2 * n, 2), np.arange(1, 2 * n, 2)])
segments = pv.PolyData(ends, lines=lines)
segments.point_data["radius"] = np.repeat(0.5 * diameters, 2)
segments.point_data["angle to z (deg)"] = np.repeat(
    np.degrees(np.arccos(np.clip(np.abs(axes[:, 2]), 0.0, 1.0))), 2
)
tubes = segments.tube(scalars="radius", absolute=True, n_sides=8)
tubes.set_active_scalars("angle to z (deg)")

pl = pv.Plotter()
pl.add_mesh(tubes, cmap="Oranges", smooth_shading=True)
pl.add_mesh(pv.Box(bounds=segments.bounds), style="wireframe", color="#24445C")
pl.show()
