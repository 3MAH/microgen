from numpy import cos, sin
import numpy as np

def gyroid(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
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


def schwarzP(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
    .. math::
       cos(x) + cos(y) + cos(z) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.schwarzP,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return cos(x) + cos(y) + cos(z)


def schwarzD(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
    .. math::
       sin(x) sin(y) sin(z) + sin(x) cos(y) cos(z) + cos(x) sin(y) cos(z) + cos(x) cos(y) sin(z) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.schwarzD,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return (
        sin(x) * sin(y) * sin(z) +
        sin(x) * cos(y) * cos(z) +
        cos(x) * sin(y) * cos(z) +
        cos(x) * cos(y) * sin(z)
    )


def neovius(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
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
    return (
        3 * cos(x) + cos(y) + cos(z) +
        4 * cos(x) * cos(y) * cos(z)
    )


def schoenIWP(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
    .. math::
       2 (cos(x) cos(y) + cos(y) cos(z) + cos(z) cos(x)) - (cos(2 x) + cos(2 y) + cos(2 z)) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.schoenIWP,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return (
        2 * (cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x)) -
        (cos(2 * x) + cos(2 * y) + cos(2 * z))
    )


def schoenFRD(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
    .. math::
       4 cos(x) cos(y) cos(z) - (cos(2 x) cos(2 y) + cos(2 y) cos(2 z) + cos(2 z) cos(2 x)) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.schoenFRD,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return (
        4 * cos(x) * cos(y) * cos(z) -
        (cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x))
    )


def fischerKochS(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
    .. math::
       cos(2 x) sin(y) cos(z) + cos(x) cos(2 y) sin(z) + sin(x) cos(y) cos(2 z) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.fischerKochS,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return (
        cos(2 * x) * sin(y) * cos(z) +
        cos(x) * cos(2 * y) * sin(z) +
        sin(x) * cos(y) * cos(2 * z)
    )


def pmy(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
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
        2 * cos(x) * cos(y) * cos(z) +
        sin(2 * x) * sin(y) +
        sin(x) * sin(2 * z) +
        sin(2 * y) * sin(z)
    )


def honeycomb(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
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
    """
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
    """
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

def honeycomb_gyroid(x: float, y: float, z: float) -> float:
    """
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


def honeycomb_schwarzP(x: float, y: float, z: float) -> float:
    """
    .. math::
       cos(x) + cos(y) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.honeycomb_schwarzP,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return cos(x) + cos(y)

def honeycomb_schwarzD(x: float, y: float, z: float) -> float:
    """
    .. math::
       cos(x) cos(y) + sin(x) sin(y) + sin(x) cos(y) + cos(x) sin(y) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.honeycomb_schwarzD,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return cos(x) * cos(y) + sin(x) * sin(y) + sin(x) * cos(y) + cos(x) * sin(y)

def honeycomb_schoenIWP(x: float, y: float, z: float) -> float:
    """
    .. math::
       cos(x) cos(y) + cos(y) + cos(x) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.surface_functions.honeycomb_schoenIWP,
           offset=0.3,
           resolution=30,
       )
       shape = geometry.sheet

       shape.plot(color='white')
    """
    return cos(x) * cos(y) + cos(y) + cos(x)

def honeycomb_lidinoid(x: float, y: float, z: float) -> float:
    """
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
    return (
        1.1 * (sin(2 * x) * cos(y) + sin(2 * y) * sin(x) + cos(x) * sin(y))
        - (cos(2 * x) * cos(2 * y) + cos(2 * y) + cos(2 * x))
    )
