"""
Generic triply periodic surfaces
================================

Any set of integer Fourier modes gives a triply periodic surface, minimal or
not.  Back row: a gyroid, and the same gyroid perturbed by a second harmonic.
Front row: a cubic surface built from two shells of modes, and a surface with
random amplitudes and phases.  All of them go through :class:`~microgen.Tpms`
like the closed-form surface functions.
"""

from itertools import product

import numpy as np
import pyvista as pv

from microgen import Tpms
from microgen.shape.surface_functions import fourier_series

gyroid_modes = [(1, 1, 0), (1, -1, 0), (0, 1, 1), (0, 1, -1), (1, 0, 1), (-1, 0, 1)]
gyroid_amplitudes = [0.5] * 6
gyroid_phases = [-0.5 * np.pi] * 6

gyroid = fourier_series(gyroid_modes, gyroid_amplitudes, gyroid_phases)

perturbed_gyroid = fourier_series(
    modes=[*gyroid_modes, (2, 0, 0), (0, 2, 0), (0, 0, 2)],
    amplitudes=[*gyroid_amplitudes, 0.25, 0.25, 0.25],
    phases=[*gyroid_phases, 0.0, 0.0, 0.0],
)

# One mode of each (m, -m) pair is enough: both give the same cosine.
shell_110 = [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]
shell_200 = [(2, 0, 0), (0, 2, 0), (0, 0, 2)]
cubic = fourier_series(
    modes=[*shell_110, *shell_200],
    amplitudes=[1.0] * len(shell_110) + [-0.6] * len(shell_200),
)

rng = np.random.default_rng(1)
low_modes = [m for m in product(range(-1, 2), repeat=3) if m > (0, 0, 0)]
random_surface = fourier_series(
    modes=low_modes,
    amplitudes=rng.normal(size=len(low_modes)),
    phases=rng.uniform(0.0, 2.0 * np.pi, size=len(low_modes)),
)

meshes = []
for i, surface in enumerate([gyroid, perturbed_gyroid, cubic, random_surface]):
    elem = Tpms(surface_function=surface, density=0.3, resolution=50)
    mesh = elem.sheet
    mesh.translate([1.2 * (i // 2), 1.2 * (i % 2), 0], inplace=True)
    meshes.append(mesh)

surfaces_block = pv.MultiBlock(meshes)

# %%
surfaces_block.plot(smooth_shading=True, color="#ef7d2d")
