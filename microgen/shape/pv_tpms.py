"""
=============================================
TPMS (:mod:`microgen.shape.tpms`)
=============================================

.. jupyter-execute::
   :hide-code:

   pyvista.global_theme.smooth_shading = True

"""
import logging
from typing import Callable, List, Union

import cadquery as cq
import numpy as np
import pyvista as pv
from numpy import cos, sin
from tqdm import tqdm

from microgen.shape.basicGeometry import BasicGeometry

SHEET = "sheet"
SKELETAL = "skeletal"
ALL = "all"

CYLINDRICAL = "cylindrical"

Field = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]

class Tpms(BasicGeometry):
    """
    Class to generate Triply Periodical Minimal Surfaces (TPMS)
    geometry from a given mathematical function, with given thickness

    functions available :
        - :class:`~microgen.shape.tpms.gyroid`
        - :class:`~microgen.shape.tpms.schwarzP`
        - :class:`~microgen.shape.tpms.schwarzD`
        - :class:`~microgen.shape.tpms.neovius`
        - :class:`~microgen.shape.tpms.schoenIWP`
        - :class:`~microgen.shape.tpms.schoenFRD`
        - :class:`~microgen.shape.tpms.fischerKochS`
        - :class:`~microgen.shape.tpms.pmy`
        - :class:`~microgen.shape.tpms.honeycomb`
        - :class:`~microgen.shape.tpms.lidinoid`
        - :class:`~microgen.shape.tpms.split_p`
    """

    def __init__(
        self,
        surface_function: Field,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        offset: Union[float, Field] = 0.,
        cell_size: Union[float, tuple[float], np.ndarray] = np.array([1., 1., 1.]),
        repeat_cell: Union[int, tuple[int], np.ndarray] = np.array([1, 1, 1]),
        resolution: int = 20,
        coordinate_system: str = "cartesian",
        cylinder_radius: float = 1.
    ) -> None:
        """
        Class used to generate TPMS geometries (sheet or skeletals parts).
        TPMS are created by default in a cube.
        The geometry of the cube can be modified using 'cell_size' parameter.
        The number of repetitions in each direction of the created geometry can be modified with 'repeat_cell' parameter.

        It is also possible to create a cylindrical geometry using 'cylindrical' coordinate system.
        For cylindrical TPMS, directions of repeat_cell must be taken as $\left(\rho, \theta, z\right)

        :param center: center of the geometry
        :param orientation: orientation of the geometry
        :param surface_function: tpms function or custom function (f(x, y, z, t) = 0)
        :param type_part: 'sheet', 'skeletal' or 'all'
        :param offset: offset of the isosurface to generate thickness
        :param cell_size: float or list of float for each dimension to set unit cell dimensions
        :param repeat_cell: integer or list of integers to repeat the geometry in each dimension
        :param resolution: unit cell resolution of the grid to compute tpms scalar fields
        :param coordinate_system: if set to 'cylindrical' then creates a cylindrical tpms, else a cartesian one is created
        :param cylinder_radius: if coordinate_system set to 'cylindrical', use cylinder_radius to define the radius of the cylinder on which the TPMS will be created        
        """
        super().__init__(shape="TPMS", center=center, orientation=orientation)

        self.surface_function = surface_function
        self.offset = offset

        self.grid = pv.StructuredGrid()
        self._sheet = None
        self._upper_skeletal = None
        self._lower_skeletal = None
        self._surface = None

        if isinstance(cell_size, (float, int)):
            self.cell_size = np.array([cell_size, cell_size, cell_size])
        else:
            self.cell_size = np.array(cell_size)

        if isinstance(repeat_cell, int):
            self.repeat_cell = np.array([repeat_cell, repeat_cell, repeat_cell])
        else:
            self.repeat_cell = np.array(repeat_cell)

        self.resolution = resolution

        self.coordinate_system = coordinate_system
        if coordinate_system == CYLINDRICAL:
            self.cylinder_radius = cylinder_radius

            unit_theta = self.cell_size[1] / cylinder_radius
            n_repeat_to_full_circle = int(2 * np.pi / unit_theta)
            self.unit_theta = 2 * np.pi / n_repeat_to_full_circle
            if self.repeat_cell[1] > n_repeat_to_full_circle:
                logging.warning(f"{n_repeat_to_full_circle} cells repeated in circular direction")
                self.repeat_cell[1] = n_repeat_to_full_circle

    @property
    def sheet(self) -> pv.PolyData:
        """
        Returns sheet part
        """
        if self._sheet is not None:
            return self._sheet

        if self.grid.dimensions == (0, 0, 0):
            self._compute_tpms_field()

        return (self.grid.clip_scalar(scalars="upper_surface", invert=False)
                         .clip_scalar(scalars="lower_surface")
                         .extract_surface())

    @property
    def upper_skeletal(self) -> pv.PolyData:
        """
        Returns upper skeletal part
        """
        if self._upper_skeletal is not None:
            return self._upper_skeletal

        if self.grid.dimensions == (0, 0, 0):
            self._compute_tpms_field()

        return self.grid.clip_scalar(scalars="upper_surface").extract_surface()

    @property
    def lower_skeletal(self) -> pv.PolyData:
        """
        Returns lower skeletal part
        """
        if self._upper_skeletal is not None:
            return self._upper_skeletal

        if self.grid.dimensions == (0, 0, 0):
            self._compute_tpms_field()

        return self.grid.clip_scalar(scalars="lower_surface", invert=False).extract_surface()

    @property
    def skeletals(self) -> tuple[pv.PolyData, pv.PolyData]:
        """
        Returns both skeletal parts
        """
        return (self.upper_skeletal, self.lower_skeletal)

    @property
    def surface(self) -> tuple[pv.PolyData]:
        """
        Returns isosurface f(x, y, z) = 0
        """
        if self._surface is not None:
            return self._surface

        if self.grid.dimensions == (0, 0, 0):
            self._compute_tpms_field()

        mesh = self.grid.contour(isosurfaces=[0.], scalars="surface")
        mesh.smooth(n_iter=100, inplace=True)
        mesh.clean(inplace=True)
        return mesh


    def _compute_tpms_field(self):
        linspaces: List[np.ndarray] = [
            np.linspace(
                -0.5 * cell_size_axis * repeat_cell_axis,
                0.5 * cell_size_axis * repeat_cell_axis,
                self.resolution * repeat_cell_axis,
            )
            for repeat_cell_axis, cell_size_axis in zip(self.repeat_cell, self.cell_size)
]

        x, y, z = np.meshgrid(*linspaces)

        k_x, k_y, k_z = 2. * np.pi / self.cell_size

        surface_function = self.surface_function(k_x * x, k_y * y, k_z * z)

        offset: Union[float, Field] = self.offset
        if not isinstance(self.offset, float):
            offset = self.offset(x, y, z)

        if self.coordinate_system == CYLINDRICAL:
            rho = x + self.cylinder_radius
            theta = y * self.unit_theta

            x = rho * np.cos(theta)
            y = rho * np.sin(theta)

        self.grid = pv.StructuredGrid(x, y, z)

        self.grid["surface"] = surface_function.ravel(order="F")
        self.grid["lower_surface"] = (surface_function - 0.5 * offset).ravel(order="F")
        self.grid["upper_surface"] = (surface_function + 0.5 * offset).ravel(order="F")

    def _create_shell(
        self,
        mesh: pv.PolyData,
        verbose: bool
    ) -> cq.Shell:
        list_of_triangles = mesh.faces.reshape(-1, 4)[:, 1:]
        list_of_triangles = np.c_[list_of_triangles, list_of_triangles[:, 0]]

        faces = []
        for i in tqdm(range(len(list_of_triangles)), disable=(not verbose)):
            ixs = list_of_triangles[i]
            lines = []
            for v1, v2 in zip(ixs[:], ixs[1:]):
                vertice_coords1 = mesh.points[v1]
                vertice_coords2 = mesh.points[v2]
                lines.append(
                    cq.Edge.makeLine(
                        cq.Vector(*vertice_coords1), cq.Vector(*vertice_coords2)
                    )
                )
            wire = cq.Wire.assembleEdges(lines)
            faces.append(cq.Face.makeFromWires(wire))

        if verbose:
            logging.warning("\nGenerating shell surface\n")
        return cq.Shell.makeShell(faces)

    def generate(self, type_part="sheet", verbose=True) -> cq.Shell:
        """
        :param type_part: part of the TPMS desired ('sheet', 'lower skeletal' or 'upper skeletal')
        :param verbose: display progressbar of the conversion to CadQuery object

        :return: CadQuery Shell object of the required TPMS part
        """
        if type_part == "sheet":
            return self._create_shell(mesh=self.sheet, verbose=verbose)
        if type_part == "lower skeletal":
            return self._create_shell(mesh=self.lower_skeletal, verbose=verbose)
        if type_part == "upper skeletal":
            return self._create_shell(mesh=self.upper_skeletal, verbose=verbose)

        logging.error(
            f"'type_part' ({type_part}) must be 'sheet', 'lower skeletal' or 'upper skeletal'"
        )
        return cq.Shell(None)


def gyroid(x: float, y: float, z: float) -> float:
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
    return (
        sin(x) * cos(y)
        + sin(y) * cos(z)
        + sin(z) * cos(x)
    )


def schwarzP(x: float, y: float, z: float) -> float:
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


def schwarzD(x: float, y: float, z: float) -> float:
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
    a = sin(x) * sin(y) * sin(z)
    b = sin(x) * cos(y) * cos(z)
    c = cos(x) * sin(y) * cos(z)
    d = cos(x) * cos(y) * sin(z)
    return a + b + c + d


def neovius(x: float, y: float, z: float) -> float:
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
    a = 3 * cos(x) + cos(y) + cos(z)
    b = 4 * cos(x) * cos(y) * cos(z)

    return a + b


def schoenIWP(x: float, y: float, z: float) -> float:
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
    a = 2 * (
        cos(x) * cos(y)
        + cos(y) * cos(z)
        + cos(z) * cos(x)
    )
    b = cos(2 * x) + cos(2 * y) + cos(2 * z)

    return a - b


def schoenFRD(x: float, y: float, z: float) -> float:
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
    a = 4 * cos(x) * cos(y) * cos(z)
    b = (
        cos(2 * x) * cos(2 * y)
        + cos(2 * y) * cos(2 * z)
        + cos(2 * z) * cos(2 * x)
    )
    return a - b


def fischerKochS(x: float, y: float, z: float) -> float:
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
    a = cos(2 * x) * sin(y) * cos(z)
    b = cos(x) * cos(2 * y) * sin(z)
    c = sin(x) * cos(y) * cos(2 * z)

    return a + b + c


def pmy(x: float, y: float, z: float) -> float:
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
    a = 2 * cos(x) * cos(y) * cos(z)
    b = sin(2 * x) * sin(y)
    c = sin(x) * sin(2 * z)
    d = sin(2 * y) * sin(z)

    return a + b + c + d


def honeycomb(x: float, y: float, z: float) -> float:
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


def lidinoid(x: float, y: float, z: float) -> float:
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
    return 0.5 * (sin(2 * x) * cos(y) * sin(z) +
                  sin(2 * y) * cos(z) * sin(x) +
                  sin(2 * z) * cos(x) * sin(y)) - \
           0.5 * (cos(2 * x) * cos(2 * y) +
                  cos(2 * y) * cos(2 * z) +
                  cos(2 * z) * cos(2 * x)) + 0.3


def split_p(x: float, y: float, z: float) -> float:
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
    return 1.1 * (sin(2 * x) * cos(y) * sin(z) +
                  sin(2 * y) * cos(z) * sin(x) +
                  sin(2 * z) * cos(x) * sin(y)) - \
           0.2 * (cos(2 * x) * cos(2 * y) +
                  cos(2 * y) * cos(2 * z) +
                  cos(2 * z) * cos(2 * x)) - \
           0.4 * (cos(2 * x) + cos(2 * y) + cos(2 * z))
