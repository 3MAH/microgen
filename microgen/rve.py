"""Representative Volume Element (RVE).

Frozen, immutable container for the RVE bounding box and its periodicity flags.

``Rve.box`` is a :class:`~microgen.cad.CadShape` wrapping an OCCT box; it is
built lazily on first access and requires the ``[cad]`` extra.  ``Rve.grid``
returns a :class:`pyvista.StructuredGrid` aligned with the cell — used by
implicit shapes to sample SDF/level-set fields.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

_DIM = 3

if TYPE_CHECKING:
    from collections.abc import Sequence

    import numpy.typing as npt
    import pyvista as pv

    from .cad import CadShape

    Vector3DType = (
        tuple[float, float, float] | Sequence[float] | npt.NDArray[np.float64]
    )
    ResolutionType = int | tuple[int, int, int] | Sequence[int]


def _validate_center(center: object) -> np.ndarray:
    if isinstance(center, (tuple, list)) and len(center) == _DIM:
        return np.asarray(center, dtype=float)
    if isinstance(center, np.ndarray) and center.shape == (_DIM,):
        return center.astype(float, copy=True)
    err_msg = f"center must be an array or Sequence of length {_DIM}"
    raise ValueError(err_msg)


def _validate_dim(dim: object) -> np.ndarray:
    if isinstance(dim, (int, float)) and not isinstance(dim, bool):
        arr = np.array([dim] * _DIM, dtype=float)
    elif isinstance(dim, (tuple, list)) and len(dim) == _DIM:
        arr = np.asarray(dim, dtype=float)
    elif isinstance(dim, np.ndarray) and dim.shape == (_DIM,):
        arr = dim.astype(float, copy=True)
    else:
        err_msg = f"dim must be an array or Sequence of length {_DIM}"
        raise ValueError(err_msg)
    if np.any(arr <= 0):
        err_msg = f"dimensions of the RVE must be greater than 0, got {arr.tolist()}"
        raise ValueError(err_msg)
    return arr


def _validate_pbc(pbc: object) -> tuple[bool, bool, bool]:
    if isinstance(pbc, bool):
        return (pbc, pbc, pbc)
    if isinstance(pbc, (tuple, list)) and len(pbc) == _DIM:
        return (bool(pbc[0]), bool(pbc[1]), bool(pbc[2]))
    err_msg = f"pbc must be a bool or Sequence[bool] of length {_DIM}"
    raise ValueError(err_msg)


@dataclass(frozen=True, init=False, eq=False, repr=False)
class Rve:
    """Representative Volume Element (RVE) — frozen.

    :param center: center of the RVE
    :param dim: dimensions of the RVE (scalar → cube)
    :param pbc: periodic-boundary-condition flags per axis ``(x, y, z)``;
        defaults to fully periodic. A single ``bool`` is broadcast to all axes.
    """

    center: np.ndarray
    dim: np.ndarray
    pbc: tuple[bool, bool, bool]

    def __init__(
        self: Rve,
        center: Vector3DType = (0, 0, 0),
        dim: float | Vector3DType = 1,
        pbc: bool | tuple[bool, bool, bool] | Sequence[bool] = (True, True, True),
    ) -> None:
        """Initialize the RVE."""
        object.__setattr__(self, "center", _validate_center(center))
        object.__setattr__(self, "dim", _validate_dim(dim))
        object.__setattr__(self, "pbc", _validate_pbc(pbc))

    @cached_property
    def min_point(self: Rve) -> np.ndarray:
        """Min corner ``center - 0.5 * dim``."""
        return self.center - 0.5 * self.dim

    @cached_property
    def max_point(self: Rve) -> np.ndarray:
        """Max corner ``center + 0.5 * dim``."""
        return self.center + 0.5 * self.dim

    @cached_property
    def box(self: Rve) -> CadShape:
        """Return a :class:`~microgen.cad.CadShape` box of the RVE (cached).

        Requires the ``[cad]`` extra.
        """
        from .cad import make_box  # noqa: PLC0415

        return make_box(tuple(self.dim), tuple(self.center))

    def grid(self: Rve, resolution: ResolutionType) -> pv.StructuredGrid:
        """Return a structured grid aligned with the RVE.

        :param resolution: points per axis — int (broadcast to all 3 axes) or
            length-3 sequence.  The grid spans ``min_point`` → ``max_point``
            with endpoints included on every axis.
        """
        import pyvista as pv  # noqa: PLC0415

        if isinstance(resolution, int):
            nx, ny, nz = resolution, resolution, resolution
        elif isinstance(resolution, (tuple, list)) and len(resolution) == _DIM:
            nx, ny, nz = (int(r) for r in resolution)
        elif isinstance(resolution, np.ndarray) and resolution.shape == (_DIM,):
            nx, ny, nz = (int(r) for r in resolution.tolist())
        else:
            err_msg = f"resolution must be an int or Sequence of length {_DIM}"
            raise ValueError(err_msg)

        xi = np.linspace(self.min_point[0], self.max_point[0], nx)
        yi = np.linspace(self.min_point[1], self.max_point[1], ny)
        zi = np.linspace(self.min_point[2], self.max_point[2], nz)
        x, y, z = np.meshgrid(xi, yi, zi, indexing="ij")
        return pv.StructuredGrid(x, y, z)

    def __repr__(self: Rve) -> str:
        return (
            f"Rve(center={tuple(self.center.tolist())}, "
            f"dim={tuple(self.dim.tolist())}, pbc={self.pbc})"
        )

    @classmethod
    def from_min_max(
        cls: type[Rve],
        min_point: Vector3DType = (-0.5, -0.5, -0.5),
        max_point: Vector3DType = (0.5, 0.5, 0.5),
        pbc: bool | tuple[bool, bool, bool] | Sequence[bool] = (True, True, True),
    ) -> Rve:
        """Generate a Rve from min and max corner points.

        :param min_point: ``(x_min, y_min, z_min)`` corner of the RVE
        :param max_point: ``(x_max, y_max, z_max)`` corner of the RVE
        :param pbc: periodic-boundary-condition flags (see ``__init__``)
        """
        lo = np.asarray(min_point, dtype=float)
        hi = np.asarray(max_point, dtype=float)
        return cls(
            center=tuple(0.5 * (lo + hi)),
            dim=tuple(np.abs(hi - lo)),
            pbc=pbc,
        )
