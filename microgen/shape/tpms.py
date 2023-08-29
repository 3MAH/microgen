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

from microgen.shape.basicGeometry import BasicGeometry
from ..operations import fuseShapes, rescale, repeatShape
from ..rve import Rve

logging.basicConfig(level=logging.INFO)
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
            self.grid.clip_scalar(scalars="upper_surface")
            .clip_scalar(scalars="lower_surface", invert=False).clean()
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
            scalars="upper_surface", invert=False
        ).clean()
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
            scalars="lower_surface"
        ).clean()
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
        return mesh

    def _create_grid(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> pv.StructuredGrid:
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

        if isinstance(self.offset, float):
            self.offset = self.offset
        elif isinstance(self.offset, Callable):
            self.offset = self.offset(x, y, z)

        self.grid["surface"] = surface_function.ravel(order="F")
        self.grid["lower_surface"] = (surface_function + 0.5 * self.offset).ravel(order="F")
        self.grid["upper_surface"] = (surface_function - 0.5 * self.offset).ravel(order="F")

    def _create_shell(self, mesh: pv.PolyData, verbose: bool) -> cq.Shell:
        triangles = mesh.faces.reshape(-1, 4)[:, 1:]
        triangles = np.c_[triangles, triangles[:, 0]]

        faces = []
        for tri in triangles:
            lines = [
                cq.Edge.makeLine(
                    cq.Vector(*mesh.points[start]), cq.Vector(*mesh.points[end])
                )
                for start, end in zip(tri[:], tri[1:])
            ]
            wire = cq.Wire.assembleEdges(lines)
            faces.append(cq.Face.makeFromWires(wire))

        if verbose:
            logging.info("\nGenerating shell surface\n")
        return cq.Shell.makeShell(faces)

    def _create_surface(
        self,
        isovalue: float = 0,
        smoothing: int = 0,
        verbose: bool = False,
    ) -> cq.Shell:
        if self.grid.dimensions == (0, 0, 0):
            self._compute_tpms_field()

        mesh = self.grid.contour(isosurfaces=[isovalue], scalars="surface")
        mesh.smooth(n_iter=smoothing, feature_smoothing=True, inplace=True)
        mesh.clean(inplace=True)

        try:
            return self._create_shell(mesh=mesh, verbose=verbose)
        except ValueError as exc:
            logging.error("Cannot create shell, try to use a higher smoothing value: %s", exc)

    def _create_surfaces(
        self,
        isovalues: list[float],
        smoothing: int = 0,
        verbose: bool = False
    ) -> list[cq.Shell]:
        """
        Create TPMS surfaces for the corresponding isovalue, return a list of cq.Shell

        :param isovalues: list of isovalues corresponding to the required surfaces
        :param smoothing: smoothing loop iterations
        :param verbose: display progressbar of the conversion to CadQuery object
        """
        shells = []
        for i, isovalue in enumerate(isovalues):
            if verbose:
                logging.info(
                    "\nGenerating surface (%d/%d) for \
                        isovalue %.2f\n",
                        i + 1, len(isovalues), isovalue
                )
            shell = self._create_surface(isovalue=isovalue, smoothing=smoothing, verbose=verbose)
            shells.append(shell)

        return shells

    def generate(
        self,
        type_part: Literal["sheet", "lower skeletal", "upper skeletal", "surface"] = "sheet",
        smoothing: int = 0,
        verbose: bool=True,
    ) -> cq.Shape:
        """
        :param type_part: part of the TPMS desired \
            ('sheet', 'lower skeletal', 'upper skeletal' or 'surface')
        :param smoothing: smoothing loop iterations
        :param verbose: display progressbar of the conversion to CadQuery object

        :return: CadQuery Shape object of the required TPMS part
        """
        if type_part not in ["sheet", "lower skeletal", "upper skeletal", "surface"]:
            raise ValueError(
                f"'type_part' ({type_part}) must be 'sheet', \
                    'lower skeletal', 'upper skeletal' or 'surface'",
            )

        if type_part == "sheet" and self.offset == 0.0:
            raise ValueError(
                "offset must be greater than 0 to generate 'sheet' part"
            )

        if type_part == "surface":
            logging.warning("offset is ignored for 'surface' part")
            return self._create_surface(isovalue=0, smoothing=smoothing, verbose=verbose)

        if not isinstance(self.offset, float):
            raise NotImplementedError(
                "Graded offset is not supported yet with the `generate` function"
            )

        eps = self.offset / 3.0
        if self.offset == 0.0:
            eps = 0.1 * np.min(self.cell_size)

        box = cq.Workplane("front").box(1, 1, 1)
        if type_part == "sheet":
            isovalues = [
                -self.offset / 2.0,
                -self.offset / 2.0 + eps,
                self.offset / 2.0 - eps,
                self.offset / 2.0,
            ]

            shells = self._create_surfaces(isovalues, smoothing, verbose)

            lower_surface = shells[0]
            lower_test_surface = shells[1]
            upper_test_surface = shells[2]
            upper_surface = shells[3]

            surface = lower_surface.fuse(upper_surface)
            test_surface = lower_test_surface.fuse(upper_test_surface)

        elif type_part == "lower skeletal":
            isovalues = [
                -self.offset / 2.0,
                -self.offset / 2.0 - eps,
            ]

            shells = self._create_surfaces(isovalues, smoothing, verbose)

            surface = shells[0]
            test_surface = shells[1]

        elif type_part == "upper skeletal":
            isovalues = [
                self.offset / 2.0,
                self.offset / 2.0 + eps,
            ]

            shells = self._create_surfaces(isovalues, smoothing, verbose)

            surface = shells[0]
            test_surface = shells[1]

        splitted_box = box.split(surface)
        tpms_solids = splitted_box.solids().all()

        # split each solid with the test surface to identify to what part type the solid belongs to
        list_solids = [
            (solid.split(test_surface).solids().size(), solid.val())
            for solid in tpms_solids
        ]

        # if the number of shapes is greater than 1, it means that the solid is splitted
        # so it belongs to the required part
        part_solids = [solid for (number, solid) in list_solids if number > 1]
        part_shapes = [cq.Shape(solid.wrapped) for solid in part_solids]
        shape = fuseShapes(cqShapeList=part_shapes, retain_edges=False) # True or False ?

        if not np.array_equal(self.cell_size, np.array([1.0, 1.0, 1.0])):
            shape = rescale(shape=shape, scale=self.cell_size)

        if not np.array_equal(self.repeat_cell, np.array([1, 1, 1])):
            shape = repeatShape(
                unit_geom=shape,
                rve=Rve(
                    dim_x=self.cell_size[0],
                    dim_y=self.cell_size[1],
                    dim_z=self.cell_size[2],
                    center=self.center,
                ),
                grid=self.repeat_cell,
            )

        return shape

    def generateVtk(
        self,
        type_part: Literal["sheet", "lower skeletal", "upper skeletal", "surface"]
    ) -> pv.PolyData:
        """
        :param type_part: part of the TPMS desireds
        :return: VTK PolyData object of the required TPMS part
        """
        if type_part == "sheet":
            return self.sheet.triangulate()
        if type_part == "lower skeletal":
            return self.lower_skeletal.triangulate()
        if type_part == "upper skeletal":
            return self.upper_skeletal.triangulate()
        if type_part == "surface":
            return self.surface.triangulate()
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
        self.cell_size[1] = self.unit_theta * radius
        if self.repeat_cell[1] == 0 or self.repeat_cell[1] > n_repeat_to_full_circle:
            logging.info(
                "%d cells repeated in circular direction", n_repeat_to_full_circle
            )
            self.repeat_cell[1] = n_repeat_to_full_circle

    def _create_grid(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> pv.StructuredGrid:
        rho = x + self.cylinder_radius
        theta = y * self.unit_theta

        return pv.StructuredGrid(rho * np.cos(theta), rho * np.sin(theta), z)


class SphericalTpms(Tpms):
    """
    Class used to generate spherical TPMS geometries (sheet or skeletals parts).
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
        Directions of cell_size and repeat_cell must be taken as the spherical coordinate system $\left(r, \theta, \phi\right)$.

        The $\theta$ and $\phi$ components of cell_size are automatically updated to the closest values that matches the spherical periodicity of the TPMS.
        If the $\theta$ or $\phi$ components of repeat_cell are 0 or greater than the periodicity of the TPMS, they are automatically set the correct number to make the full sphere.

        :param radius: radius of the sphere on which the center of the TPMS is located
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

        self.sphere_radius = radius

        unit_theta = self.cell_size[1] / radius
        n_repeat_theta_to_join = int(np.pi / unit_theta)
        self.unit_theta = np.pi / n_repeat_theta_to_join
        self.cell_size[1] = self.unit_theta * radius # true only on theta = pi/2
        if self.repeat_cell[1] == 0 or self.repeat_cell[1] > n_repeat_theta_to_join:
            logging.info(
                "%d cells repeated in theta direction", n_repeat_theta_to_join
            )
            self.repeat_cell[1] = n_repeat_theta_to_join

        unit_phi = self.cell_size[2] / radius
        n_repeat_phi_to_join = int(2 * np.pi / unit_phi)
        self.unit_phi = 2 * np.pi / n_repeat_phi_to_join
        self.cell_size[2] = self.unit_phi * radius
        if self.repeat_cell[2] == 0 or self.repeat_cell[2] > n_repeat_phi_to_join:
            logging.info(
                "%d cells repeated in phi direction", n_repeat_phi_to_join
            )
            self.repeat_cell[2] = n_repeat_phi_to_join

    def _create_grid(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> pv.StructuredGrid:
        rho = x + self.sphere_radius
        theta = y * self.unit_theta + np.pi / 2.0
        phi = z * self.unit_phi

        return pv.StructuredGrid(
            rho * np.sin(theta) * np.cos(phi),
            rho * np.sin(theta) * np.sin(phi),
            rho * np.cos(theta),
        )
