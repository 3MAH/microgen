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
           surface_function=microgen.tpms.gyroid,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

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
           surface_function=microgen.tpms.schwarzP,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

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
           surface_function=microgen.tpms.schwarzD,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

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
           surface_function=microgen.tpms.neovius,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

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
           surface_function=microgen.tpms.schoenIWP,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

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
           surface_function=microgen.tpms.schoenFRD,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

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
           surface_function=microgen.tpms.fischerKochS,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

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
           surface_function=microgen.tpms.pmy,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

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
           surface_function=microgen.tpms.honeycomb,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white')
    """
    return sin(x) * cos(y) + sin(y) + cos(z)


def lidinoid(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """
    .. math::
       0.5 * (sin(2 * x) * cos(y) * sin(z) +
              sin(2 * y) * cos(z) * sin(x) +
              sin(2 * z) * cos(x) * sin(y)) -
       0.5 * (cos(2 * x) * cos(2 * y) +
              cos(2 * y) * cos(2 * z) +
              cos(2 * z) * cos(2 * x)) + 0.3 = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.lidinoid,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

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
       1.1 * (sin(2 * x) * cos(y) * sin(z) +
              sin(2 * y) * cos(z) * sin(x) +
              sin(2 * z) * cos(x) * sin(y)) -
       0.2 * (cos(2 * x) * cos(2 * y) +
              cos(2 * y) * cos(2 * z) +
              cos(2 * z) * cos(2 * x)) -
       0.4 * (cos(2 * x) + cos(2 * y) + cos(2 * z)) = 0

    .. jupyter-execute::
       :hide-code:

       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.split_p,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

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
