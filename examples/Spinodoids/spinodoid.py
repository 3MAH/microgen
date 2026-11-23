"""
Spinodoid
=========

Spinodoids mimic the bicontinuous microstructures of spinodal
decomposition. A :class:`microgen.Spinodoid` is a Gaussian random field: an
FFT of white noise filtered by a Gaussian shell of radius ``k0`` (the
dominant wave number) and width ``bandwidth``. The solid is the part where
the field is above a threshold; giving ``density`` makes microgen find the
threshold that matches this solid volume fraction.

Each kept Fourier mode lives on the reciprocal lattice of the cell, so the
field is exactly periodic: the spinodoid is a ready-made periodic RVE.
"""

# sphinx_gallery_thumbnail_number = 1
import numpy as np
import pyvista as pv

from microgen import Phase, Spinodoid

# %%
# An isotropic spinodoid with 50 % of solid. ``seed`` makes the random field
# reproducible.
spinodoid = Spinodoid(k0=30.0, bandwidth=4.0, density=0.5, resolution=64, seed=1)
volume_fraction = spinodoid.grid_solid.volume / np.prod(spinodoid.cell_size)
print(f"Solid volume fraction: {volume_fraction:.3f}")

# %%
# The field is periodic: its values on opposite faces of the cell match.
y, z = np.meshgrid(np.linspace(0, 1, 20), np.linspace(0, 1, 20))
left = spinodoid.evaluate(np.zeros(y.size), y.ravel(), z.ravel())
right = spinodoid.evaluate(np.ones(y.size), y.ravel(), z.ravel())
print(f"Max gap between the x = 0 and x = 1 faces: {np.abs(left - right).max():.2e}")

# %%
# A field-backed :class:`microgen.Phase` gives the connected components and
# the moments without CAD. At 50 % both phases percolate: nearly all the
# solid is one connected piece.
phase = Phase.from_shape(spinodoid, resolution=64)
volumes = sorted((piece.volume for piece in phase.pieces), reverse=True)
print(f"Pieces: {len(volumes)}, largest holds {volumes[0] / sum(volumes):.1%}")
print(f"Centre of mass: {phase.center_of_mass}")

surface = spinodoid.generate_surface_mesh()
surface.save("spinodoid.stl")
surface.plot(color="#ef7d2d", smooth_shading=True)

# %%
# ``anisotropy`` stretches the wave vectors along each axis before filtering.
# A large factor along an axis removes the modes with a component along it:
# ``(1, 1, 3)`` gives columns along z, ``(3, 3, 1)`` gives lamellae normal
# to z.
variants = {
    "isotropic": spinodoid,
    "columnar (1, 1, 3)": Spinodoid(
        k0=30.0, bandwidth=4.0, density=0.5, anisotropy=(1, 1, 3), resolution=64, seed=1
    ),
    "lamellar (3, 3, 1)": Spinodoid(
        k0=30.0, bandwidth=4.0, density=0.5, anisotropy=(3, 3, 1), resolution=64, seed=1
    ),
}

pl = pv.Plotter(shape=(1, 3), window_size=(1500, 500))
for i, (name, shape) in enumerate(variants.items()):
    pl.subplot(0, i)
    pl.add_text(name, font_size=12)
    pl.add_mesh(shape.generate_surface_mesh(), color="#ef7d2d", smooth_shading=True)
pl.link_views()
pl.show()
