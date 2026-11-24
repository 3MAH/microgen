"""TPMS surface functions."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import autograd.numpy as np
import numpy as onp
from autograd.numpy import cos, sin

from ._frep_grf import _FrepGRF


def gyroid(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Gyroid.

    .. math::
       sin(x) cos(y) + sin(y) cos(z) + sin(z) cos(x) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.gyroid,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x)


def schwarz_p(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Schwarz P.

    .. math::
       cos(x) + cos(y) + cos(z) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.schwarz_p,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return cos(x) + cos(y) + cos(z)


def schwarz_d(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Schwarz D.

    .. math::
       sin(x) sin(y) sin(z) + sin(x) cos(y) cos(z) +\
          cos(x) sin(y) cos(z) + cos(x) cos(y) sin(z) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.schwarz_d,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return (
        sin(x) * sin(y) * sin(z)
        + sin(x) * cos(y) * cos(z)
        + cos(x) * sin(y) * cos(z)
        + cos(x) * cos(y) * sin(z)
    )


def neovius(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Neovius.

    .. math::
        3 cos(x) + cos(y) + cos(z) + 4 cos(x) cos(y) cos(z) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.neovius,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return 3 * cos(x) + cos(y) + cos(z) + 4 * cos(x) * cos(y) * cos(z)


def schoen_iwp(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Schoen IWP.

    .. math::
       2 (cos(x) cos(y) + cos(y) cos(z) + cos(z) cos(x)) \
        - (cos(2 x) + cos(2 y) + cos(2 z)) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.schoen_iwp,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return 2 * (cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x)) - (
        cos(2 * x) + cos(2 * y) + cos(2 * z)
    )


def schoen_frd(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Schoen FRD.

    .. math::
       4 cos(x) cos(y) cos(z) \
        - (cos(2 x) cos(2 y) + cos(2 y) cos(2 z) + cos(2 z) cos(2 x)) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.schoen_frd,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return 4 * cos(x) * cos(y) * cos(z) - (
        cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x)
    )


def fischer_koch_s(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Fischer-Koch S.

    .. math::
       cos(2 x) sin(y) cos(z) + cos(x) cos(2 y) sin(z) + sin(x) cos(y) cos(2 z) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.fischer_koch_s,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return (
        cos(2 * x) * sin(y) * cos(z)
        + cos(x) * cos(2 * y) * sin(z)
        + sin(x) * cos(y) * cos(2 * z)
    )


def pmy(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """PMY.

    .. math::
       2 cos(x) cos(y) cos(z) + sin(2 x) sin(y) + sin(x) sin(2 z) + sin(2 y) sin(z) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.pmy,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return (
        2 * cos(x) * cos(y) * cos(z)
        + sin(2 * x) * sin(y)
        + sin(x) * sin(2 * z)
        + sin(2 * y) * sin(z)
    )


def honeycomb(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Honneycomb.

    .. math::
       sin(x) cos(y) + sin(y) + cos(z) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.honeycomb,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return sin(x) * cos(y) + sin(y) + cos(z)


def lidinoid(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Lidinoid.

    .. math::
       0.5 (sin(2 x) cos(y) sin(z) +
            sin(2 y) cos(z) sin(x) +
            sin(2 z) cos(x) sin(y)) -
       0.5 (cos(2 x) cos(2 y) +
            cos(2 y) cos(2 z) +
            cos(2 z) cos(2 x)) + 0.3 = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.lidinoid,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return (
        0.5
        * (
            sin(2 * x) * cos(y) * sin(z)
            + sin(2 * y) * cos(z) * sin(x)
            + sin(2 * z) * cos(x) * sin(y)
        )
        - 0.5
        * (cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x))
        + 0.3
    )


def split_p(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Split P.

    .. math::
       1.1 (sin(2 x) cos(y) sin(z) +
            sin(2 y) cos(z) sin(x) +
            sin(2 z) cos(x) sin(y)) -
       0.2 (cos(2 x) cos(2 y) +
            cos(2 y) cos(2 z) +
            cos(2 z) cos(2 x)) -
       0.4 (cos(2 x) + cos(2 y) + cos(2 z)) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.split_p,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return (
        1.1
        * (
            sin(2 * x) * cos(y) * sin(z)
            + sin(2 * y) * cos(z) * sin(x)
            + sin(2 * z) * cos(x) * sin(y)
        )
        - 0.2
        * (cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x))
        - 0.4 * (cos(2 * x) + cos(2 * y) + cos(2 * z))
    )


def honeycomb_gyroid(x: np.ndarray, y: np.ndarray, _: np.ndarray) -> np.ndarray:
    """Honeycomb Gyroid.

    .. math::
       sin(x) cos(y) + sin(y) + cos(x) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.honeycomb_gyroid,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return sin(x) * cos(y) + sin(y) + cos(x)


def honeycomb_schwarz_p(x: np.ndarray, y: np.ndarray, _: np.ndarray) -> np.ndarray:
    """Honeycomb Schwarz P.

    .. math::
       cos(x) + cos(y) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.honeycomb_schwarz_p,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return cos(x) + cos(y)


def honeycomb_schwarz_d(x: np.ndarray, y: np.ndarray, _: np.ndarray) -> np.ndarray:
    """Honneycomb Schwarz D.

    .. math::
       cos(x) cos(y) + sin(x) sin(y) + sin(x) cos(y) + cos(x) sin(y) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.honeycomb_schwarz_d,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return cos(x) * cos(y) + sin(x) * sin(y) + sin(x) * cos(y) + cos(x) * sin(y)


def honeycomb_schoen_iwp(x: np.ndarray, y: np.ndarray, _: np.ndarray) -> np.ndarray:
    """Honneycomb Schoen IWP.

    .. math::
       cos(x) cos(y) + cos(y) + cos(x) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.honeycomb_schoen_iwp,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return cos(x) * cos(y) + cos(y) + cos(x)


def honeycomb_lidinoid(x: np.ndarray, y: np.ndarray, _: np.ndarray) -> np.ndarray:
    """Honeycomb Lidinoid.

    .. math::
       1.1 (sin(2 x) cos(y) + sin(2 y) sin(x) + cos(x) sin(y))
       - (cos(2 x) cos(2 y) + cos(2 y) + cos(2 x)) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.honeycomb_lidinoid,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return 1.1 * (sin(2 * x) * cos(y) + sin(2 * y) * sin(x) + cos(x) * sin(y)) - (
        cos(2 * x) * cos(2 * y) + cos(2 * y) + cos(2 * x)
    )


def fourier_series(
    modes: Sequence[Sequence[int]],
    amplitudes: Sequence[float],
    phases: Sequence[float] | None = None,
    constant: float = 0.0,
) -> Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]:
    r"""Triply periodic surface function defined by a truncated Fourier series.

    .. math::
       \sum_i A_i \cos(h_i x + k_i y + l_i z + \varphi_i) + c = 0

    Every mode :math:`(h_i, k_i, l_i)` is a triplet of integers, so the field
    is :math:`2\pi`-periodic along each axis and the returned function can be
    used as the ``surface_function`` of :class:`~microgen.shape.tpms.Tpms`,
    like the closed-form functions of this module.  The surface is triply
    periodic but, in general, not minimal.

    A sine term is a cosine of phase :math:`-\pi/2`, and a product of
    harmonics expands into sum and difference modes, e.g.
    :math:`\sin(x) \cos(y) = \tfrac{1}{2} \sin(x + y) + \tfrac{1}{2} \sin(x - y)`.

    :param modes: ``(K, 3)`` integer reciprocal-lattice indices
        :math:`(h_i, k_i, l_i)`
    :param amplitudes: ``(K,)`` amplitudes :math:`A_i`
    :param phases: ``(K,)`` phases :math:`\varphi_i` in radians (default: 0)
    :param constant: constant term :math:`c`
    :return: surface function ``f(x, y, z)``
    :raises ValueError: if a mode is not an integer triplet, or if the number
        of amplitudes or phases differs from the number of modes

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.fourier_series(
               modes=[(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (2, 1, 0)],
               amplitudes=[1.0, 1.0, 1.0, 0.8, 0.5],
               phases=[0.0, 0.0, 0.0, 1.0, -1.5],
           ),
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    modes_arr = onp.asarray(modes, dtype=float)
    if modes_arr.ndim != 2 or modes_arr.shape[1] != 3:  # noqa: PLR2004
        err_msg = f"modes must have shape (K, 3), got {modes_arr.shape}"
        raise ValueError(err_msg)
    if not onp.array_equal(modes_arr, onp.round(modes_arr)):
        err_msg = "modes must be integers for the field to be periodic"
        raise ValueError(err_msg)
    amplitudes_arr = onp.asarray(amplitudes, dtype=float)
    phases_arr = (
        onp.zeros(len(modes_arr))
        if phases is None
        else onp.asarray(phases, dtype=float)
    )
    for name, arr in (("amplitudes", amplitudes_arr), ("phases", phases_arr)):
        if arr.shape != (len(modes_arr),):
            err_msg = f"{name} must have shape ({len(modes_arr)},), got {arr.shape}"
            raise ValueError(err_msg)

    # Tpms feeds the surface function with 2*pi/cell_size-scaled coordinates,
    # hence the 2*pi period here.
    return _FrepGRF(
        modes=modes_arr.astype(onp.int64),
        amplitudes=amplitudes_arr,
        phases=phases_arr,
        threshold=-float(constant),
        cell_size=onp.full(3, 2.0 * onp.pi),
    ).evaluate
