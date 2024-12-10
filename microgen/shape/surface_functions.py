"""TPMS surface functions."""

import numpy as np
from numpy import cos, sin


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


def honeycomb_gyroid(x: float, y: float, _: float) -> float:
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


def honeycomb_schwarz_p(x: float, y: float, _: float) -> float:
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


def honeycomb_schwarz_d(x: float, y: float, _: float) -> float:
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


def honeycomb_schoen_iwp(x: float, y: float, _: float) -> float:
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


def honeycomb_lidinoid(x: float, y: float, _: float) -> float:
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


# Deprecated functions
schwarzP = schwarz_p  # noqa: N816
schwarzD = schwarz_d  # noqa: N816
schoenIWP = schoen_iwp  # noqa: N816
schoenFRD = schoen_frd  # noqa: N816
fischerKochS = fischer_koch_s  # noqa: N816
honeycomb_schwarzP = honeycomb_schwarz_p  # noqa: N816
honeycomb_schoenIWP = honeycomb_schoen_iwp  # noqa: N816
honeycomb_schwarzD = honeycomb_schwarz_d  # noqa: N816
