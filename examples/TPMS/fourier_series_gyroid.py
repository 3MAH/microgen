"""
Gyroid from its Fourier modes
=============================

The closed-form TPMS functions are short truncated Fourier series.
:func:`~microgen.shape.surface_functions.fourier_series` builds a surface
function from a list of integer modes, amplitudes and phases; here the six
modes of the gyroid give back the same sheet as
:func:`~microgen.shape.surface_functions.gyroid`.
"""

import numpy as np

from microgen import Tpms
from microgen.shape.surface_functions import fourier_series, gyroid

# sin(x) cos(y) = 1/2 sin(x + y) + 1/2 sin(x - y), and sin(u) = cos(u - pi/2)
modes = [
    (1, 1, 0),
    (1, -1, 0),
    (0, 1, 1),
    (0, 1, -1),
    (1, 0, 1),
    (-1, 0, 1),
]
gyroid_series = fourier_series(
    modes=modes,
    amplitudes=[0.5] * len(modes),
    phases=[-0.5 * np.pi] * len(modes),
)

points = np.random.default_rng(0).uniform(-np.pi, np.pi, (3, 1000))
print(f"max difference: {np.abs(gyroid_series(*points) - gyroid(*points)).max():.1e}")

geometry = Tpms(
    surface_function=gyroid_series,
    density=0.30,
    resolution=30,
)
shape = geometry.generate_surface_mesh(type_part="sheet")

# %%
shape.plot(smooth_shading=True, color="#ef7d2d")
