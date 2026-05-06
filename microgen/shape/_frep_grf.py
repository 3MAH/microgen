"""F-rep Gaussian Random Field backend.

==========================================================
F-rep GRF backend (:mod:`microgen.shape._frep_grf`)
==========================================================

Internal-only backend that stores a sparse periodic Fourier-mode expansion
and exposes it as a continuously evaluable F-rep callable.
:meth:`_FrepGRF.from_fft_grf` keeps the Fourier coefficients of a
Gaussian-Random-Field generated via FFT of white noise filtered by a
Gaussian shell — the formulation used by most spinodoid mechanics
literature, and the basis on which any future GSTools-style covariance
preset would slot in.

Periodicity is automatic: every kept mode lives on the reciprocal lattice
(integer index in each axis), so ``evaluate(x + L) == evaluate(x)``
bit-for-bit.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

_DIM = 3


def _normalize_cell_size(
    cell_size: Sequence[float] | float,
) -> npt.NDArray[np.float64]:
    """Coerce ``cell_size`` to a length-3 float ndarray."""
    if isinstance(cell_size, (int, float)):
        return np.array([cell_size, cell_size, cell_size], dtype=float)
    arr = np.asarray(cell_size, dtype=float)
    if arr.shape != (3,):
        err_msg = f"cell_size must be scalar or 3-tuple, got {arr.shape}"
        raise ValueError(err_msg)
    return arr


def _lex_positive_modes(modes: npt.NDArray[np.int64]) -> npt.NDArray[np.bool_]:
    """Return mask selecting one of each (m, -m) pair (real-field Hermitian symmetry).

    A mode is kept iff its first non-zero component is positive.
    """
    if len(modes) == 0:
        return np.zeros(0, dtype=bool)
    nonzero = modes != 0
    first_nonzero_idx = np.argmax(nonzero, axis=1)
    has_nonzero = nonzero.any(axis=1)
    rows = np.arange(len(modes))
    sign_at_first = np.where(has_nonzero, modes[rows, first_nonzero_idx], 0)
    return sign_at_first > 0


@dataclass
class _FrepGRF:
    """Sparse periodic Fourier-mode expansion of a scalar field.

    The implicit field is::

        f(x) = sum_i amplitudes[i] * cos(2*pi * modes[i] @ x / cell_size + phases[i])

    with ``modes[i]`` a 3-tuple of integers (reciprocal-lattice indices) and
    ``cell_size`` the period along each axis. The ``threshold`` is subtracted
    from ``f(x)`` so that the iso-surface ``f(x) - threshold = 0`` encloses
    the desired porosity.

    :ivar modes: ``(K, 3)`` integer array of reciprocal-lattice indices
    :ivar amplitudes: ``(K,)`` real amplitudes
    :ivar phases: ``(K,)`` phases in radians
    :ivar threshold: scalar; the iso-surface offset
    :ivar cell_size: ``(3,)`` array; period along each axis
    """

    modes: npt.NDArray[np.int64]
    amplitudes: npt.NDArray[np.float64]
    phases: npt.NDArray[np.float64]
    threshold: float
    cell_size: npt.NDArray[np.float64]

    def _field(
        self,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Evaluate the raw field ``Σᵢ Aᵢ cos(kᵢ·x + φᵢ)`` (no threshold)."""
        x_arr = np.asarray(x, dtype=np.float64)
        y_arr = np.asarray(y, dtype=np.float64)
        z_arr = np.asarray(z, dtype=np.float64)
        out_shape = np.broadcast_shapes(x_arr.shape, y_arr.shape, z_arr.shape)
        flat = np.broadcast_arrays(x_arr, y_arr, z_arr)
        points = np.stack([a.ravel() for a in flat], axis=-1)
        if self.modes.size == 0:
            return np.zeros(out_shape, dtype=np.float64)
        # Single (M, K) buffer reused for angles → cos(angles+phase) → amplitude-weighted contribs.
        scaled = points / self.cell_size
        angles = scaled @ self.modes.T
        angles *= 2.0 * np.pi
        angles += self.phases
        np.cos(angles, out=angles)
        angles *= self.amplitudes
        return angles.sum(axis=1).reshape(out_shape)

    def evaluate(
        self,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Evaluate ``f(x, y, z) - threshold`` at arbitrary points.

        Vectorized over the input arrays (any broadcastable shape).
        """
        return self._field(x, y, z) - self.threshold

    def sample_on_grid(
        self,
        resolution: int | Sequence[int],
        repeat_cell: int | Sequence[int] = 1,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Sample the field on a structured grid covering one (or more) cells.

        Uses ``linspace(0, L, N+1, endpoint=True)`` so the grid extrema land
        bit-exactly at ``L`` along each axis — required for the cap-plane
        classification in :func:`microgen.cad.mesh_to_periodic_shell`.

        :param resolution: per-axis number of *cells* (the returned grid has
            ``resolution + 1`` points per axis); int or 3-tuple.
        :param repeat_cell: tile the cell this many times along each axis.
        :return: ``(grid_points, field_values)`` where ``grid_points`` has
            shape ``(Nx+1, Ny+1, Nz+1, 3)`` and ``field_values`` has shape
            ``(Nx+1, Ny+1, Nz+1)``.
        """
        if isinstance(resolution, int):
            res = np.array([resolution, resolution, resolution], dtype=int)
        else:
            res = np.asarray(resolution, dtype=int)
            if res.shape != (3,):
                err_msg = f"resolution must be int or 3-tuple, got shape {res.shape}"
                raise ValueError(err_msg)
        if isinstance(repeat_cell, int):
            rep = np.array([repeat_cell, repeat_cell, repeat_cell], dtype=int)
        else:
            rep = np.asarray(repeat_cell, dtype=int)
        nx, ny, nz = (res * rep).tolist()
        lx, ly, lz = (self.cell_size * rep).tolist()
        xs = np.linspace(0.0, lx, nx + 1, endpoint=True)
        ys = np.linspace(0.0, ly, ny + 1, endpoint=True)
        zs = np.linspace(0.0, lz, nz + 1, endpoint=True)
        xv, yv, zv = np.meshgrid(xs, ys, zs, indexing="ij")
        values = self.evaluate(xv, yv, zv)
        points = np.stack([xv, yv, zv], axis=-1)
        return points, values

    @staticmethod
    def from_fft_grf(  # noqa: PLR0913
        k0: float,
        bandwidth: float,
        anisotropy: Sequence[float],
        cell_size: Sequence[float] | float,
        resolution: int,
        seed: int | None,
        mode_amplitude_threshold: float = 1e-3,
    ) -> _FrepGRF:
        """Build an _FrepGRF whose modes come from FFT of filtered white noise.

        Generates a Gaussian Random Field by filtering white noise with a
        Gaussian shell ``Γ(K) = exp(-(K_eff − k₀)² / (2·bandwidth²))`` in
        Fourier space, then keeps the surviving Fourier coefficients
        (``|F̂| > mode_amplitude_threshold · max(|F̂|)``) as the sparse
        F-rep mode set. ``threshold`` is left at 0; the caller sets it from
        a target porosity.

        :param k0: dominant wave number (Gaussian shell center)
        :param bandwidth: Gaussian shell width (was ``theta``)
        :param anisotropy: ``(aₓ, a_y, a_z)`` axis stretch factors applied to
            wavenumbers in the shell filter ``K_eff² = (aₓ Kₓ)² + …``
        :param cell_size: period along each axis (scalar or 3-tuple)
        :param resolution: FFT grid resolution (single int, applied per axis)
        :param seed: RNG seed; ``None`` for non-deterministic behavior
        :param mode_amplitude_threshold: keep modes with
            ``|F̂| > threshold * max(|F̂|)``. Smaller → more modes, more
            accurate F-rep but slower evaluation.
        """
        cs = _normalize_cell_size(cell_size)
        ax, ay, az = (float(a) for a in anisotropy)

        rng = np.random.default_rng(seed)
        n = int(resolution)
        idx = np.fft.fftfreq(n, d=1.0 / n).astype(np.int64)
        kx_idx, ky_idx, kz_idx = np.meshgrid(idx, idx, idx, indexing="ij")
        kx_phys = 2.0 * np.pi * kx_idx / cs[0]
        ky_phys = 2.0 * np.pi * ky_idx / cs[1]
        kz_phys = 2.0 * np.pi * kz_idx / cs[2]
        k_eff = np.sqrt((ax * kx_phys) ** 2 + (ay * ky_phys) ** 2 + (az * kz_phys) ** 2)
        gamma = np.exp(-((k_eff - k0) ** 2) / (2.0 * bandwidth**2))

        white_noise = rng.standard_normal((n, n, n))
        spectrum = np.fft.fftn(white_noise) * gamma
        mag = np.abs(spectrum)
        max_mag = float(mag.max()) or 1.0
        keep = mag > mode_amplitude_threshold * max_mag
        keep[0, 0, 0] = False  # DC term only shifts the mean
        kx_keep = kx_idx[keep]
        ky_keep = ky_idx[keep]
        kz_keep = kz_idx[keep]
        modes = np.stack([kx_keep, ky_keep, kz_keep], axis=-1).astype(np.int64)
        coeffs = spectrum[keep]
        amplitudes = np.abs(coeffs).astype(np.float64) / (n**_DIM)
        phases = np.angle(coeffs).astype(np.float64)

        # Real-field Hermitian symmetry: each (m, -m) pair contributes the same cosine,
        # so keep one and double the amplitude.
        keep_lex = _lex_positive_modes(modes)
        modes = modes[keep_lex]
        amplitudes = 2.0 * amplitudes[keep_lex]
        phases = phases[keep_lex]

        return _FrepGRF(
            modes=modes,
            amplitudes=amplitudes,
            phases=phases,
            threshold=0.0,
            cell_size=cs,
        )


def compute_threshold_for_porosity(
    grf: _FrepGRF,
    porosity: float,
    n_samples: int = 10_000,
    seed: int | None = None,
) -> float:
    """Find the iso-value ``T`` such that ``Pr(f(x) >= T) ≈ porosity``.

    Monte Carlo over the unit cell. The smooth periodic field gives a stable
    ``(1−porosity)``-quantile estimate at modest sample counts; 10k samples
    keeps the threshold within ~0.5 % of the asymptotic value for typical
    spinodoid spectra.
    """
    rng = np.random.default_rng(seed)
    pts = rng.uniform(0.0, 1.0, size=(n_samples, 3)) * grf.cell_size
    values = grf._field(pts[:, 0], pts[:, 1], pts[:, 2])
    return float(np.quantile(values, 1.0 - porosity))
