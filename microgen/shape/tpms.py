"""
=============================================
TPMS (:mod:`microgen.shape.tpms`)
=============================================

.. jupyter-execute::
   :hide-code:

   pyvista.global_theme.smooth_shading = True

"""
import logging
from typing import Callable, List, Union, Sequence, Literal

import cadquery as cq
import numpy as np
import pyvista as pv
# from tqdm import tqdm

from microgen.shape.basicGeometry import BasicGeometry
from ..operations import fuseShapes, rescale, repeatShape
# from ..rve import Rve

Field = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]


class Tpms(BasicGeometry):
    """
    Class to generate Triply Periodical Minimal Surfaces (TPMS)
    geometry from a given mathematical function, with given offset

    functions available :
        - :class:`~microgen.shape.surface_functions.gyroid`
        - :class:`~microgen.shape.surface_functions.schwarzP`
        - :class:`~microgen.shape.surface_functions.schwarzD`
        - :class:`~microgen.shape.surface_functions.neovius`
        - :class:`~microgen.shape.surface_functions.schoenIWP`
        - :class:`~microgen.shape.surface_functions.schoenFRD`
        - :class:`~microgen.shape.surface_functions.fischerKochS`
        - :class:`~microgen.shape.surface_functions.pmy`
        - :class:`~microgen.shape.surface_functions.honeycomb`
        - :class:`~microgen.shape.surface_functions.lidinoid`
        - :class:`~microgen.shape.surface_functions.split_p`
    """

    def __init__(
        self,
        surface_function: Field,
        offset: Union[float, Field] = 0.0,
        phase_shift: Sequence[float] = (0.0, 0.0, 0.0),
        cell_size: Union[float, Sequence[float]] = 1.0,
        repeat_cell: Union[int, Sequence[int]] = 1,
        resolution: int = 20,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
    ) -> None:
        """
        Class used to generate TPMS geometries (sheet or skeletals parts).
        TPMS are created by default in a cube.
        The geometry of the cube can be modified using 'cell_size' parameter.
        The number of repetitions in each direction of the created geometry can be modified with the 'repeat_cell' parameter.

        :param center: center of the geometry
        :param orientation: orientation of the geometry
        :param surface_function: tpms function or custom function (f(x, y, z) = 0)
        :param offset: offset of the isosurface to generate thickness
        :param phase_shift: phase shift of the isosurface $f(x + \phi_x, y + \phi_y, z + \phi_z, t) = 0$
        :param cell_size: float or list of float for each dimension to set unit cell dimensions
        :param repeat_cell: integer or list of integers to repeat the geometry in each dimension
        :param resolution: unit cell resolution of the grid to compute tpms scalar fields
        """
        super().__init__(shape="TPMS", center=center, orientation=orientation)

        self.surface_function = surface_function
        self.offset = offset
        self.phase_shift = phase_shift

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

        self.grid = self._create_grid(x, y, z)

        k_x, k_y, k_z = 2.0 * np.pi / self.cell_size
        surface_function = self.surface_function(
            k_x * (x + self.phase_shift[0]),
            k_y * (y + self.phase_shift[1]),
            k_z * (z + self.phase_shift[2])
        )

        offset: Union[float, np.ndarray] = 0.0
        if isinstance(self.offset, float):
            offset = self.offset
        elif isinstance(self.offset, Callable):
            offset = self.offset(x, y, z)

        self.grid["surface"] = surface_function.ravel(order="F")
        self.grid["lower_surface"] = (surface_function - 0.5 * offset).ravel(order="F")
        self.grid["upper_surface"] = (surface_function + 0.5 * offset).ravel(order="F")

    def _create_shell(self, mesh: pv.PolyData, verbose: bool) -> cq.Shell:
        triangles = mesh.faces.reshape(-1, 4)[:, 1:]
        triangles = np.c_[triangles, triangles[:, 0]]

        faces = []
        # for i in tqdm(range(len(triangles)), disable=not verbose):
        for i in range(len(triangles)):
            tri = triangles[i]
            lines = [
                cq.Edge.makeLine(
                    cq.Vector(*mesh.points[start]), cq.Vector(*mesh.points[end])
                )
                for start, end in zip(tri[:], tri[1:])
            ]
            wire = cq.Wire.assembleEdges(lines)
            faces.append(cq.Face.makeFromWires(wire))

        if verbose:
            logging.warning("\nGenerating shell surface\n")
        return cq.Shell.makeShell(faces)
    
    def _create_surfaces(
        self,
        isovalues: list[float] = [0],
    ) -> list[cq.Shell]:
        """
        Create TPMS surfaces for the corresponding isovalue, return a list of cq.Shell

        :param numsber_surfaces: number of surfaces
        :param isovalues: height isovalues of the given tpms function
        :param nSample: surface file name
        :param smoothing: smoothing loop iterations
        """
        self._compute_tpms_field()

        shells = []
        for isovalue in isovalues:
            mesh = self.grid.contour(isosurfaces=[isovalue], scalars="surface")
            mesh.smooth(n_iter=100, inplace=True)
            mesh.clean(inplace=True)

            shell = self._create_shell(mesh=mesh, verbose=False)
            shells.append(shell)

        return shells

    def generate(self, type_part="sheet", verbose=True) -> cq.Shape:
        """
        :param type_part: part of the TPMS desired ('sheet', 'lower skeletal', 'upper skeletal' or 'surface')
        :param verbose: display progressbar of the conversion to CadQuery object

        :return: CadQuery Shape object of the required TPMS part
        """
        if type_part not in ["sheet", "lower skeletal", "upper skeletal", "surface"]:
            raise ValueError(
                f"'type_part' ({type_part}) must be 'sheet', \
                    'lower skeletal', 'upper skeletal' or 'surface'",
            )

        if type_part == "surface":
            shell = self._create_shell(mesh=self.surface, verbose=verbose)
            return cq.Shape(shell.wrapped)

        isovalues = [
            -self.offset,
            -self.offset / 3.0,
            self.offset / 3.0,
            self.offset,
        ]
        shells = self._create_surfaces(isovalues=isovalues)

        face_cut_tp = shells[2]
        face_cut_tm = shells[1]
        face_cut_p = shells[3]
        face_cut_m = shells[0]

        box_wp = cq.Workplane("front").box(1, 1, 1)

        box_cut_wp = box_wp.split(face_cut_p)
        box_cut_wp = box_cut_wp.split(face_cut_m)

        box_workplanes = box_cut_wp.solids().all()  # type: list[cq.Workplane]

        list_shapes = [
            (wp.split(face_cut_tp).split(face_cut_tm).solids().size(), wp.val())
            for wp in box_workplanes
        ]
        if type_part == "sheet":
            sheet = [shape for (number, shape) in list_shapes if number > 1]
            to_fuse = [cq.Shape(shape.wrapped) for shape in sheet]
            shape = fuseShapes(to_fuse, True)
        elif type_part == "skeletal":
            skeletal = [shape for (number, shape) in list_shapes if number == 1]
            to_fuse = [cq.Shape(shape.wrapped) for shape in skeletal]
            shape = fuseShapes(to_fuse, False)

        return shape
        # raise NotImplementedError("Only 'surface' is implemented for now")

    def generateVtk(
        self,
        type_part: Literal["sheet", "lower skeletal", "upper skeletal", "surface"]
    ) -> pv.PolyData:
        """
        :param type_part: part of the TPMS desired
        :return: VTK PolyData object of the required TPMS part
        """
        if type_part == "sheet":
            return self.sheet
        if type_part == "lower skeletal":
            return self.lower_skeletal
        if type_part == "upper skeletal":
            return self.upper_skeletal
        if type_part == "surface":
            return self.surface

        raise ValueError(
            f"type_part ({type_part}) must be 'sheet', \
                'lower skeletal', 'upper skeletal' or 'surface'"
        )


class CylindricalTpms(Tpms):
    """
    Class used to generate cylindrical TPMS geometries (sheet or skeletals parts).
    """
    def __init__(
        self,
        radius: float,
        surface_function: Field,
        offset: Union[float, Field] = 0.0,
        phase_shift: Sequence[float] = (0.0, 0.0, 0.0),
        cell_size: Union[float, Sequence[float]] = 1.0,
        repeat_cell: Union[int, Sequence[int]] = 1,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        resolution: int = 20,
    ):
        """
        Directions of cell_size and repeat_cell must be taken as the cylindrical coordinate system $\left(\rho, \theta, z\right)$.

        The $\theta$ component of cell_size is automatically updated to the closest value that matches the cylindrical periodicity of the TPMS.
        If the $\theta$ component of repeat_cell is 0 or greater than the periodicity of the TPMS, it is automatically set the correct number to make the full cylinder.

        :param radius: radius of the cylinder on which the center of the TPMS is located
        :param surface_function: tpms function or custom function (f(x, y, z) = 0)
        :param offset: offset of the isosurface to generate thickness
        :param phase_shift: phase shift of the tpms function $f(x + \phi_x, y + \phi_y, z + \phi_z) = 0$
        :param cell_size: float or list of float for each dimension to set unit cell dimensions
        :param repeat_cell: integer or list of integers to repeat the geometry in each dimension
        :param center: center of the geometry
        :param orientation: orientation of the geometry
        :param resolution: unit cell resolution of the grid to compute tpms scalar fields
        """
        super().__init__(
            surface_function=surface_function,
            offset=offset,
            phase_shift=phase_shift,
            cell_size=cell_size,
            repeat_cell=repeat_cell,
            resolution=resolution,
            center=center,
            orientation=orientation,
        )

        self.cylinder_radius = radius

        unit_theta = self.cell_size[1] / radius
        n_repeat_to_full_circle = int(2 * np.pi / unit_theta)
        self.unit_theta = 2 * np.pi / n_repeat_to_full_circle
        self.cell_size[1] = self.unit_theta * radius # TODO : check if this is correct
        if self.repeat_cell[1] == 0 or self.repeat_cell[1] > n_repeat_to_full_circle:
            logging.warning(
                "%d cells repeated in circular direction", n_repeat_to_full_circle
            )
            self.repeat_cell[1] = n_repeat_to_full_circle

    def _create_grid(self, x, y, z):
        rho = x + self.cylinder_radius
        theta = y * self.unit_theta

        return pv.StructuredGrid(rho * np.cos(theta), rho * np.sin(theta), z)
