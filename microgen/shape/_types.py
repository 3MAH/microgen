"""Shared type aliases for the shape package.

Centralising these here avoids duplicate definitions across ``shape.py``,
``tpms.py`` and downstream modules.  All implicit shapes use the same
``Field`` callable and ``BoundsType`` AABB representation.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import numpy.typing as npt

# Implicit scalar field: ``(x, y, z) -> array``, with negative values inside.
Field = Callable[
    [npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]],
    npt.NDArray[np.float64],
]

# Axis-aligned bounding box: ``(xmin, xmax, ymin, ymax, zmin, zmax)``.
BoundsType = tuple[float, float, float, float, float, float]

# Period of an intrinsically-periodic shape: ``(Lx, Ly, Lz)``.
# When set, ``field(x + Lx, y, z) == field(x, y, z)`` (and analogously for y, z).
# A ``None`` period means the shape is not intrinsically periodic.
PeriodType = tuple[float, float, float]

__all__ = ["BoundsType", "Field", "PeriodType"]
