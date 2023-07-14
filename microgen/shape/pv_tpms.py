"""
=============================================
TPMS (:mod:`microgen.shape.tpms`)
=============================================

.. jupyter-execute::
   :hide-code:

   pyvista.global_theme.smooth_shading = True

"""
import logging
from typing import Callable, List, Union, Sequence

import cadquery as cq
import numpy as np
import pyvista as pv
from numpy import cos, sin
from tqdm import tqdm

from microgen.shape.basicGeometry import BasicGeometry

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
        offset: Union[float, Field] = 0.0,
        cell_size: Union[float, Sequence[float]] = 1.0,
        repeat_cell: Union[int, Sequence[int]] = 1,
        resolution: int = 20,
    ) -> None:
        """
        Class used to generate TPMS geometries (sheet or skeletals parts).
        TPMS are created by default in a cube.
        The geometry of the cube can be modified using 'cell_size' parameter.
        The number of repetitions in each direction of the created geometry can be modified with the 'repeat_cell' parameter.

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
        """
        super().__init__(shape="TPMS", center=center, orientation=orientation)

        self.surface_function = surface_function
        self.offset = offset

        self.grid = pv.StructuredGrid()
        self._sheet: pv.PolyData = None
        self._upper_skeletal: pv.PolyData = None
        self._lower_skeletal: pv.PolyData = None
        self._surface: pv.PolyData = None

        if isinstance(cell_size, (float, int)):
            self.cell_size = np.array([cell_size, cell_size, cell_size])
        elif len(cell_size) == 3:
            self.cell_size = np.array(cell_size)
        else:
            raise ValueError("cell_size must be a float or a sequence of 3 floats")

        if isinstance(repeat_cell, int):
            self.repeat_cell = np.array([repeat_cell, repeat_cell, repeat_cell])
        elif len(repeat_cell) == 3:
            self.repeat_cell = np.array(repeat_cell)
        else:
            raise ValueError("repeat_cell must be an int or a sequence of 3 ints")

        self.resolution = resolution

    @property
    def sheet(self) -> pv.PolyData:
        """
        Returns sheet part
        """
        if self._sheet is not None:
            return self._sheet

        if self.grid.dimensions == (0, 0, 0):
            self._compute_tpms_field()

        self._sheet: pv.PolyData = (
            self.grid.clip_scalar(scalars="upper_surface", invert=False)
            .clip_scalar(scalars="lower_surface")
            .extract_surface()
        )
        return self._sheet

    @property
    def upper_skeletal(self) -> pv.PolyData:
        """
        Returns upper skeletal part
        """
        if self._upper_skeletal is not None:
            return self._upper_skeletal

        if self.grid.dimensions == (0, 0, 0):
            self._compute_tpms_field()

        self._upper_skeletal: pv.PolyData = self.grid.clip_scalar(
            scalars="upper_surface"
        ).extract_surface()
        return self._upper_skeletal

    @property
    def lower_skeletal(self) -> pv.PolyData:
        """
        Returns lower skeletal part
        """
        if self._lower_skeletal is not None:
            return self._lower_skeletal

        if self.grid.dimensions == (0, 0, 0):
            self._compute_tpms_field()

        self._lower_skeletal: pv.PolyData = self.grid.clip_scalar(
            scalars="lower_surface", invert=False
        ).extract_surface()
        return self._lower_skeletal

    @property
    def skeletals(self) -> tuple[pv.PolyData, pv.PolyData]:
        """
        Returns both skeletal parts
        """
        return (self.upper_skeletal, self.lower_skeletal)

    @property
    def surface(self) -> pv.PolyData:
        """
        Returns isosurface f(x, y, z) = 0
        """
        if self._surface is not None:
            return self._surface

        if self.grid.dimensions == (0, 0, 0):
            self._compute_tpms_field()

        mesh: pv.PolyData = self.grid.contour(isosurfaces=[0.0], scalars="surface")
        mesh.smooth(n_iter=100, inplace=True)
        mesh.clean(inplace=True)
        return mesh

    def _create_grid(self, x, y, z):
        return pv.StructuredGrid(x, y, z)

    def _compute_tpms_field(self):
        linspaces: List[np.ndarray] = [
            np.linspace(
                -0.5 * cell_size_axis * repeat_cell_axis,
                0.5 * cell_size_axis * repeat_cell_axis,
                self.resolution * repeat_cell_axis,
            )
            for repeat_cell_axis, cell_size_axis in zip(
                self.repeat_cell, self.cell_size
            )
        ]

        x, y, z = np.meshgrid(*linspaces)

        k_x, k_y, k_z = 2.0 * np.pi / self.cell_size

        self.grid = self._create_grid(x, y, z)

        surface_function = self.surface_function(k_x * x, k_y * y, k_z * z)

        offset: Union[float, np.ndarray] = 0.0
        if isinstance(self.offset, float):
            offset = self.offset
        elif isinstance(self.offset, Callable):
            offset = self.offset(x, y, z)

        self.grid["surface"] = surface_function.ravel(order="F")
        self.grid["lower_surface"] = (surface_function - 0.5 * offset).ravel(order="F")
        self.grid["upper_surface"] = (surface_function + 0.5 * offset).ravel(order="F")

    def _create_shell(self, mesh: pv.PolyData, verbose: bool) -> cq.Shell:
        list_of_triangles = mesh.faces.reshape(-1, 4)[:, 1:]
        list_of_triangles = np.c_[list_of_triangles, list_of_triangles[:, 0]]

        faces = []
        for i in tqdm(range(len(list_of_triangles)), disable=not verbose):
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
            "'type_part' (%s) must be 'sheet', 'lower skeletal' or 'upper skeletal'",
            type_part,
        )
        return cq.Shell(None)
    
class CylindricalTpms(Tpms):
    def __init__(
        self,
        radius: float,
        surface_function: Field,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        offset: Union[float, Field] = 0.0,
        cell_size: Union[float, Sequence[float]] = 1.0,
        repeat_cell: Union[int, Sequence[int]] = 1,
        resolution: int = 20,
    ):
        super().__init__(
            surface_function=surface_function,
            center=center,
            orientation=orientation,
            offset=offset,
            cell_size=cell_size,
            repeat_cell=repeat_cell,
            resolution=resolution,
        )

        self.cylinder_radius = radius

        unit_theta = self.cell_size[1] / radius
        n_repeat_to_full_circle = int(2 * np.pi / unit_theta)
        self.unit_theta = 2 * np.pi / n_repeat_to_full_circle
        if self.repeat_cell[1] == 0 or self.repeat_cell[1] > n_repeat_to_full_circle:
            logging.warning(
                "%d cells repeated in circular direction", n_repeat_to_full_circle
            )
            self.repeat_cell[1] = n_repeat_to_full_circle

    def _create_grid(self, x, y, z):
        rho = x + self.cylinder_radius
        theta = y * self.unit_theta

        return pv.StructuredGrid(rho * np.cos(theta), rho * np.sin(theta), z)
