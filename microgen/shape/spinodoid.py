"""Spinodoid.

=============================================
Spinodoid (:mod:`microgen.shape.spinodoid`)
=============================================

A Spinodoid is a periodic stochastic microstructure built as a Gaussian
Random Field generated via FFT of white noise filtered by a Gaussian
shell. The Fourier coefficients of the filtered spectrum are kept as a
sparse expansion on the reciprocal lattice, so the implicit field is
continuously evaluable everywhere and participates in implicit boolean
ops on the same footing as TPMS shapes.

The implicit field is exactly periodic on ``cell_size × repeat_cell`` —
every kept mode lives on the reciprocal lattice (integer index in each
axis), so periodicity is a data-structure invariant, not a runtime flag.
"""

from __future__ import annotations

import contextlib
from collections.abc import Sequence
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import pyvista as pv

from ..operations import rotate
from ._frep_grf import _FrepGRF, _normalize_cell_size, compute_threshold_for_porosity
from .periodic_shell import mesh_to_periodic_shell, try_make_solid
from .shape import Shape

if TYPE_CHECKING:
    from microgen.cad import CadShape
    from microgen.shape import KwargsGenerateType, Vector3DType

    from ._types import BoundsType


class Spinodoid(Shape):
    """Class to generate a periodic Spinodoid microstructure.

    Construction is in two stages:

    1. **Spectral**: an FFT of filtered white noise yields a sparse set of
       reciprocal-lattice modes ``(modes, amplitudes, phases)`` whose
       magnitudes survive a relative-amplitude threshold.
    2. **Geometric**: the iso-surface ``f(x) = offset`` defines the solid
       (``f(x) ≥ offset`` is inside). Either ``offset`` or ``density`` must
       be given, mirroring :class:`microgen.Tpms`:

       - ``offset``: the iso-value directly (a scalar in field units).
       - ``density``: target solid volume fraction in (0, 1); microgen
         finds the ``offset`` matching this fraction via Monte Carlo on
         the F-rep callable.

       The structured grid and ``grid_solid`` / ``surface`` properties are
       derived by sampling the callable on a ``(resolution+1)^3`` grid
       covering the full unit cell.

    The implicit field ``f(x) = Σᵢ Aᵢ cos(kᵢ·x + φᵢ) − offset`` is
    *bit-for-bit periodic* on ``cell_size`` (every ``kᵢ`` lives on the
    reciprocal lattice). Periodicity is an invariant — there is no
    ``periodic=True`` flag.

    :param k0: dominant wave number (Gaussian shell center).
    :param bandwidth: Gaussian shell width.
    :param offset: iso-value defining the solid surface (mutually
        exclusive with ``density``).
    :param density: target solid volume fraction in (0, 1) (mutually
        exclusive with ``offset``).
    :param anisotropy: ``(aₓ, a_y, a_z)`` axis-stretch factors applied
        to wavenumbers in the shell filter ``K_eff² = (aₓ Kₓ)² + …``.
    :param cell_size: unit-cell side length(s); scalar or 3-tuple.
    :param repeat_cell: integer tiling along each axis.
    :param resolution: FFT grid resolution per cell (the structured grid
        will have ``resolution + 1`` points per axis to cover both
        endpoints of the periodic cell exactly).
    :param mode_amplitude_threshold: keep FFT modes with
        ``|F̂| > threshold * max(|F̂|)``.
    :param seed: RNG seed; ``None`` for non-deterministic behavior.
    :param center: geometry center after construction (translation applied
        in :meth:`generate_surface_mesh` / :meth:`generate_cad`).
    :param orientation: rotation angles (Euler) applied at the end.
    """

    def __init__(  # noqa: PLR0913
        self: Spinodoid,
        k0: float = 30.0,
        bandwidth: float = 5.0,
        offset: float | None = None,
        density: float | None = None,
        anisotropy: Sequence[float] = (1.0, 1.0, 1.0),
        cell_size: float | Sequence[float] = 1.0,
        repeat_cell: int | Sequence[int] = 1,
        resolution: int = 32,
        mode_amplitude_threshold: float = 1e-3,
        seed: int | None = None,
        center: Vector3DType = (0.0, 0.0, 0.0),
        orientation: Vector3DType = (0.0, 0.0, 0.0),
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the Spinodoid object — see class docstring."""
        super().__init__(center=center, orientation=orientation, **kwargs)

        if offset is not None and density is not None:
            err_msg = (
                "offset and density cannot be given at the same time. Give only one."
            )
            raise ValueError(err_msg)
        if offset is None and density is None:
            err_msg = "offset or density must be given. Give one of them."
            raise ValueError(err_msg)
        if density is not None and not 0.0 < density < 1.0:
            err_msg = f"density must be in (0, 1). Given: {density}"
            raise ValueError(err_msg)

        self.k0 = float(k0)
        self.bandwidth = float(bandwidth)
        self.anisotropy = np.asarray(anisotropy, dtype=float)
        self.offset = float(offset) if offset is not None else None
        self.density = float(density) if density is not None else None
        self.resolution = int(resolution)
        self.mode_amplitude_threshold = float(mode_amplitude_threshold)
        self.seed = seed

        self._init_cell_parameters(cell_size, repeat_cell)

        self._frep: _FrepGRF
        self._setup_frep_field()

    def _init_cell_parameters(
        self: Spinodoid,
        cell_size: float | Sequence[float],
        repeat_cell: int | Sequence[int],
    ) -> None:
        self.cell_size = _normalize_cell_size(cell_size)

        if isinstance(repeat_cell, int):
            self.repeat_cell = np.array(
                [repeat_cell, repeat_cell, repeat_cell],
                dtype=int,
            )
        elif len(repeat_cell) == 3:
            self.repeat_cell = np.asarray(repeat_cell, dtype=int)
        else:
            err_msg = f"`repeat_cell` must have length 3. Given: {repeat_cell}"
            raise ValueError(err_msg)

    def _setup_frep_field(self: Spinodoid) -> None:
        """Build the F-rep modes, resolve the iso-value, wire ``_func``/``_bounds``."""
        self._frep = _FrepGRF.from_fft_grf(
            k0=self.k0,
            bandwidth=self.bandwidth,
            anisotropy=self.anisotropy,
            cell_size=self.cell_size,
            resolution=self.resolution,
            seed=self.seed,
            mode_amplitude_threshold=self.mode_amplitude_threshold,
        )

        if self.offset is not None:
            self._frep.threshold = self.offset
        else:
            self._frep.threshold = compute_threshold_for_porosity(
                self._frep,
                porosity=self.density,
                seed=self.seed,
            )

        # Shape uses "negative = inside"; _frep.evaluate is positive inside, so flip the sign.
        frep = self._frep

        def _signed_field(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            return -frep.evaluate(x, y, z)

        self._field = _signed_field
        lx, ly, lz = (self.cell_size * self.repeat_cell).tolist()
        self._bounds = (0.0, float(lx), 0.0, float(ly), 0.0, float(lz))
        # Spinodoid's field is *bit-exact* periodic on ``cell_size`` along each
        # axis — every kept Fourier mode lives on the reciprocal lattice, so
        # ``frep.evaluate(p + cell_size) == frep.evaluate(p)`` (see _frep_grf.py).
        cs = np.asarray(self.cell_size, dtype=float)
        self._period = (float(cs[0]), float(cs[1]), float(cs[2]))

    @cached_property
    def grid(self: Spinodoid) -> pv.StructuredGrid:
        """Return the structured-grid sampling of the implicit field."""
        points, values = self._frep.sample_on_grid(
            resolution=self.resolution,
            repeat_cell=self.repeat_cell,
        )
        sg = pv.StructuredGrid(points[..., 0], points[..., 1], points[..., 2])
        sg["surface"] = values.ravel(order="F")
        return sg

    @cached_property
    def grid_solid(self: Spinodoid) -> pv.UnstructuredGrid:
        """Return the volumetric solid grid (cells where ``f >= 0``)."""
        return self.grid.clip_scalar(scalars="surface", value=0.0, invert=False)

    @cached_property
    def surface(self: Spinodoid) -> pv.PolyData:
        """Return the iso-surface mesh of the spinodoid."""
        return self.grid_solid.extract_surface().clean().triangulate()

    # ---- public generators -------------------------------------------------

    def generate_volume_mesh(
        self: Spinodoid,
        **_: KwargsGenerateType,
    ) -> pv.UnstructuredGrid:
        """Generate the volumetric solid grid (3D cells) with center+rotation applied."""
        grid_vol = self.grid_solid.copy()
        grid_vol = rotate(grid_vol, center=(0, 0, 0), rotation=self.orientation)
        return grid_vol.translate(xyz=self.center)

    def generate_surface_mesh(
        self: Spinodoid,
        bounds: BoundsType | None = None,
        resolution: int | None = None,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate the iso-surface as a PolyData with center+rotation applied.

        The native path uses the cached structured-grid surface whose resolution
        is fixed at construction time (``Spinodoid(resolution=…)``). When
        ``bounds`` or ``resolution`` is explicitly provided here, fall back to
        the implicit (marching-cubes-on-SDF) base implementation so polymorphic
        callers can re-sample the field at an arbitrary resolution / bbox.
        """
        if bounds is not None or resolution is not None:
            return super().generate_surface_mesh(
                bounds=bounds,
                resolution=resolution if resolution is not None else 50,
            )
        polydata = self.surface.copy()
        polydata = rotate(polydata, center=(0, 0, 0), rotation=self.orientation)
        return polydata.translate(xyz=self.center)

    def generate_cad(
        self: Spinodoid,
        **_: KwargsGenerateType,
    ) -> CadShape:
        """Generate an OCCT CAD shape (Solid or Compound) of the spinodoid.

        Mirrors :meth:`microgen.Tpms.generate_cad`: extracts the structured-grid
        surface, builds a closed multi-face periodic shell with planar BREP
        faces on the cell sides plus per-triangle planar faces on the
        spinodoid interior, and applies the ``center`` / ``orientation``
        rigid transform. Spinodoids commonly decompose into multiple disjoint
        solid components — the return type may be a Compound of fused Solids.

        Requires the optional ``[cad]`` install extra.
        """
        mesh = self.surface
        pts = np.asarray(mesh.points, dtype=np.float64)
        tris = mesh.faces.reshape(-1, 4)[:, 1:].astype(np.int64)

        shape = try_make_solid(mesh_to_periodic_shell(pts, tris, self._bounds))
        shape = rotate(obj=shape, center=(0, 0, 0), rotation=self.orientation)
        shape = shape.translate(self.center)
        # Fallback for OCCT Volume() on solids it flags invalid (rigid transforms preserve volume).
        with contextlib.suppress(AttributeError, ValueError):
            shape._mesh_volume = float(abs(self.grid_solid.volume))
        return shape
