"""
=============================================
TPMS (:mod:`microgen.shape.tpms`)
=============================================

.. jupyter-execute::
   :hide-code:

   pyvista.global_theme.smooth_shading = True
   pyvista.global_theme.split_sharp_edges = True

"""
import logging
from typing import Callable, List, Literal, Optional, Sequence, Union

import cadquery as cq
import numpy as np
import pyvista as pv
from scipy.optimize import root_scalar

from microgen.shape.basicGeometry import BasicGeometry

from ..operations import fuseShapes, repeatShape, rescale, rotateEuler, rotatePvEuler
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
        - :class:`~microgen.shape.surface_functions.honeycomb_gyroid`
        - :class:`~microgen.shape.surface_functions.honeycomb_schwarzP`
        - :class:`~microgen.shape.surface_functions.honeycomb_schwarzD`
        - :class:`~microgen.shape.surface_functions.honeycomb_schoenIWP`
        - :class:`~microgen.shape.surface_functions.honeycomb_lidinoid`
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
        density: Optional[float] = None,
    ) -> None:
        """
        Class used to generate TPMS geometries (sheet or skeletals parts).
        TPMS are created by default in a cube.
        The geometry of the cube can be modified using 'cell_size' parameter.
        The number of repetitions in each direction of the created geometry \
            can be modified with the 'repeat_cell' parameter.

        :param center: center of the geometry
        :param orientation: orientation of the geometry
        :param surface_function: tpms function or custom function (f(x, y, z) = 0)
        :param offset: offset of the isosurface to generate thickness
        :param phase_shift: phase shift of the isosurface \
            $f(x + \\phi_x, y + \\phi_y, z + \\phi_z, t) = 0$
        :param cell_size: float or list of float for each dimension to set unit cell dimensions
        :param repeat_cell: integer or list of integers to repeat the geometry in each dimension
        :param resolution: unit cell resolution of the grid to compute tpms scalar fields
        :param density: density percentage of the generated geometry (0 < density < 1) \
            If density is given, the offset is automatically computed to fit the density \
                (performance is slower than when using the offset)
        """
        super().__init__(shape="TPMS", center=center, orientation=orientation)

        self.surface_function = surface_function
        self.offset = offset
        self.phase_shift = phase_shift

        self.grid = pv.StructuredGrid()
        self._sheet = None
        self._upper_skeletal = None
        self._lower_skeletal = None
        self._surface = None

        self._init_cell_parameters(cell_size, repeat_cell)

        self.resolution = resolution
        self._compute_tpms_field()

        if density is not None and not 0.0 < density <= 1.0:
            raise ValueError("density must be between 0 and 1")
        self.density = density

    def _max_density(
        self,
        part_type: Literal["sheet", "lower skeletal", "upper skeletal"],
        resolution: Optional[int] = None,
    ) -> float:
        if part_type == "sheet":
            return 1.0
        tpms = Tpms(
            surface_function=self.surface_function,
            offset=0.0,
            resolution=resolution if resolution is not None else self.resolution,
        )
        return tpms.generateVtk(type_part=part_type).volume / tpms.grid.volume

    @classmethod
    def offset_from_density(
        cls,
        surface_function: Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray],
        part_type: Literal["sheet", "lower skeletal", "upper skeletal"],
        density: float | Literal["max"] = "max",
        resolution: int = 20,
    ) -> float:
        """
        Returns the offset corresponding to the required density for the specified part of the given surface_function.

        :param surface_function: tpms function
        :param part_type: type of the part (sheet, lower skeletal or upper skeletal)
        :param density: Required density, 0.5 for 50%
        :param resolution: resolution of the tpms used to compute the offset

        :return: corresponding offset value
        """
        if not isinstance(density, (int, float)) and density != "max":
            raise ValueError("density must be a float between 0 and 1 or 'max'")
        if density == "max":
            if part_type == "sheet":
                tpms = Tpms(surface_function=surface_function, resolution=resolution)
                return 2.0 * np.max(tpms.grid["surface"])
            return 0.0  # skeletal

        tpms = Tpms(surface_function=surface_function, density=density)
        tpms._compute_offset_to_fit_density(part_type=part_type, resolution=resolution)

        if isinstance(tpms.offset, float):
            return tpms.offset
        raise ValueError("offset must be a float")

    def _compute_offset_to_fit_density(
        self,
        part_type: Literal["sheet", "lower skeletal", "upper skeletal"],
        resolution: Optional[int] = None,
    ) -> None:
        if self.density is None:
            raise ValueError("density must be given to compute offset")
        temp_tpms = Tpms(
            surface_function=self.surface_function,
            resolution=resolution if resolution is not None else self.resolution,
        )
        max_density = temp_tpms._max_density(part_type=part_type, resolution=resolution)
        if self.density > max_density:
            raise ValueError(
                f"density ({self.density}) must be lower than {max_density} for \
                    the {part_type} part of the given TPMS function"
            )

        if self.density == max_density:
            offset = 2.0 * np.max(self.grid["surface"]) if part_type == "sheet" else 0.0
            self._update_offset(offset)
            return

        bound = 2.0 * np.max(self.grid["surface"])
        grid_volume = temp_tpms.grid.volume

        polydata_func = getattr(temp_tpms, f"vtk_{part_type.replace(' ', '_')}")

        def density(offset: float) -> float:
            temp_tpms._update_offset(offset)
            return polydata_func().volume / grid_volume

        computed_offset = root_scalar(
            lambda offset: density(offset) - self.density,
            bracket=[-bound, bound],
            method="secant",
            x0=0.5,
        ).root
        self._update_offset(computed_offset)
        logging.info("computed offset = %.3f", computed_offset)

    def _init_cell_parameters(
        self,
        cell_size: Union[float, Sequence[float]],
        repeat_cell: Union[int, Sequence[int]],
    ):
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

    def vtk_sheet(self) -> pv.PolyData:
        return self.grid.clip_scalar(scalars="upper_surface").clip_scalar(
            scalars="lower_surface", invert=False
        )

    def vtk_upper_skeletal(self) -> pv.PolyData:
        return self.grid.clip_scalar(scalars="upper_surface", invert=False)

    def vtk_lower_skeletal(self) -> pv.PolyData:
        return self.grid.clip_scalar(scalars="lower_surface")

    @property
    def sheet(self) -> pv.PolyData:
        """
        Returns sheet part
        """
        if self._sheet is not None:
            return self._sheet

        if self.density is not None:
            self._compute_offset_to_fit_density(part_type="sheet")

        self._sheet = self.vtk_sheet().clean().triangulate()
        return self._sheet

    @property
    def upper_skeletal(self) -> pv.PolyData:
        """
        Returns upper skeletal part
        """
        if self._upper_skeletal is not None:
            return self._upper_skeletal

        if self.density is not None:
            self._compute_offset_to_fit_density(part_type="upper skeletal")

        self._upper_skeletal = self.vtk_upper_skeletal().clean().triangulate()
        return self._upper_skeletal

    @property
    def lower_skeletal(self) -> pv.PolyData:
        """
        Returns lower skeletal part
        """
        if self._lower_skeletal is not None:
            return self._lower_skeletal

        if self.density is not None:
            self._compute_offset_to_fit_density(part_type="lower skeletal")

        self._lower_skeletal = self.vtk_lower_skeletal().clean().triangulate()
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

        self._surface = self.grid.contour(
            isosurfaces=[0.0], scalars="surface"
        ).triangulate()
        return self._surface

    def _create_grid(
        self, x: np.ndarray, y: np.ndarray, z: np.ndarray
    ) -> pv.StructuredGrid:
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
        tpms_field = self.surface_function(
            k_x * (x + self.phase_shift[0]),
            k_y * (y + self.phase_shift[1]),
            k_z * (z + self.phase_shift[2]),
        )

        self.grid["surface"] = tpms_field.ravel(order="F")
        self._update_offset(self.offset)

    def _update_offset(self, offset: Union[float, Callable]) -> None:
        if isinstance(offset, float):
            self.offset = offset
        elif isinstance(offset, Callable):
            self.offset = offset(self.grid.x, self.grid.y, self.grid.z).ravel("F")

        self.grid["lower_surface"] = self.grid["surface"] + 0.5 * self.offset
        self.grid["upper_surface"] = self.grid["surface"] - 0.5 * self.offset

    def _create_shell(self, mesh: pv.PolyData, verbose: bool) -> cq.Shell:
        if not mesh.is_all_triangles:
            mesh.triangulate(inplace=True)  # useless ?
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
        isovalue: Union[float, np.ndarray] = 0.0,
        smoothing: int = 0,
        verbose: bool = False,
    ) -> cq.Shell:
        if isinstance(isovalue, (int, float)):
            scalars = self.grid["surface"] - isovalue
        elif isinstance(isovalue, np.ndarray):
            scalars = self.grid["surface"] - isovalue.ravel(order="F")

        mesh = self.grid.contour(isosurfaces=[0.0], scalars=scalars)
        mesh.smooth(n_iter=smoothing, feature_smoothing=True, inplace=True)
        mesh.clean(inplace=True)

        try:
            shell = self._create_shell(mesh=mesh, verbose=verbose)
        except ValueError as exc:
            raise ValueError(
                f"Cannot create shell, try to use a higher smoothing value: {exc}"
            ) from exc
        return shell

    def _create_surfaces(
        self, isovalues: list[float], smoothing: int = 0, verbose: bool = False
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
                logging.info("\nGenerating surface (%d/%d)", i + 1, len(isovalues))
            shell = self._create_surface(
                isovalue=isovalue, smoothing=smoothing, verbose=verbose
            )
            shells.append(shell)

        return shells

    def _generate_sheet_surfaces(
        self, eps: float, smoothing: int, verbose: bool
    ) -> tuple[cq.Shape, cq.Shape]:
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
        return surface, test_surface

    def _generate_lower_skeletal_surfaces(
        self, eps: float, smoothing: int, verbose: bool
    ) -> tuple[cq.Shape, cq.Shape]:
        isovalues = [
            -self.offset / 2.0,
            -self.offset / 2.0 - eps,
        ]

        shells = self._create_surfaces(isovalues, smoothing, verbose)

        surface = shells[0]
        test_surface = shells[1]
        return surface, test_surface

    def _generate_upper_skeletal_surfaces(
        self, eps: float, smoothing: int, verbose: bool
    ) -> tuple[cq.Shape, cq.Shape]:
        isovalues = [
            self.offset / 2.0,
            self.offset / 2.0 + eps,
        ]

        shells = self._create_surfaces(isovalues, smoothing, verbose)

        surface = shells[0]
        test_surface = shells[1]
        return surface, test_surface

    def _extract_part_from_box(
        self,
        type_part: Literal["sheet", "lower skeletal", "upper skeletal", "surface"],
        eps: float,
        smoothing: int,
        verbose: bool,
    ):
        box = cq.Workplane("front").box(1, 1, 1)

        surface, test_surface = getattr(
            self, f"_generate_{type_part.replace(' ', '_')}_surfaces"
        )(eps, smoothing, verbose)

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
        return fuseShapes(
            cqShapeList=part_shapes, retain_edges=False
        )  # True or False ?

    def generate(
        self,
        type_part: Literal[
            "sheet", "lower skeletal", "upper skeletal", "surface"
        ] = "sheet",
        smoothing: int = 0,
        verbose: bool = True,
        algo_resolution: Optional[int] = None,
    ) -> cq.Shape:
        """
        :param type_part: part of the TPMS desired \
            ('sheet', 'lower skeletal', 'upper skeletal' or 'surface')
        :param smoothing: smoothing loop iterations
        :param verbose: display progressbar of the conversion to CadQuery object
        :param algo_resolution: if offset must be computed to fit density, \
            resolution of the temporary TPMS used to compute the offset

        :return: CadQuery Shape object of the required TPMS part
        """
        if type_part not in ["sheet", "lower skeletal", "upper skeletal", "surface"]:
            raise ValueError(
                f"'type_part' ({type_part}) must be 'sheet', \
                    'lower skeletal', 'upper skeletal' or 'surface'",
            )

        if type_part == "surface":
            if self.offset != 0.0:
                logging.warning("offset is ignored for 'surface' part")
            if self.density is not None:
                logging.warning("density is ignored for 'surface' part")
            return self._create_surface(
                isovalue=0, smoothing=smoothing, verbose=verbose
            )

        if self.density is not None:
            self._compute_offset_to_fit_density(
                part_type=type_part, resolution=algo_resolution
            )

        if "skeletal" in type_part:
            if (
                isinstance(self.offset, (int, float)) and self.offset < 0.0
            ):  # scalar offset = 0 is working
                raise NotImplementedError(
                    "generating 'skeletal' parts with a negative \
                        offset value is not implemented yet"
                )
            if isinstance(self.offset, np.ndarray) and np.any(self.offset <= 0.0):
                raise NotImplementedError(
                    "generating 'skeletal' parts with negative or zero \
                        offset values is not implemented yet"
                )
        elif type_part == "sheet":
            if np.any(self.offset <= 0.0):
                if np.all(self.offset <= 0.0):
                    raise ValueError(
                        "offset must be greater than 0 to generate 'sheet' part"
                    )
                raise NotImplementedError(
                    "generating 'sheet' parts with negative or zero \
                        offset values is not implemented yet"
                )

        offset_limit = 2.0 * np.max(self.grid["surface"])
        if np.all(self.offset > offset_limit):
            raise ValueError(
                f"offset ({self.offset}) must be lower \
                        than {offset_limit} for the given TPMS function"
            )

        eps = self.offset / 3.0
        if isinstance(self.offset, float) and self.offset == 0.0:
            eps = 0.1 * np.min(self.cell_size)

        shape = self._extract_part_from_box(type_part, eps, smoothing, verbose)

        density = shape.Volume() / (np.prod(self.repeat_cell) * np.prod(self.cell_size))
        logging.info(f"TPMS density = {density:.2%}")

        shape = rotateEuler(
            obj=shape,
            center=(0, 0, 0),
            psi=self.orientation[0],
            theta=self.orientation[1],
            phi=self.orientation[2],
        )
        return shape.translate(self.center)

    def generateVtk(
        self,
        type_part: Literal["sheet", "lower skeletal", "upper skeletal", "surface"],
        algo_resolution: Optional[int] = None,
    ) -> pv.PolyData:
        """
        :param type_part: part of the TPMS desireds
        :param algo_resolution: if offset must be computed to fit density, \
            resolution of the temporary TPMS used to compute the offset

        :return: VTK PolyData object of the required TPMS part
        """
        if type_part == "surface":
            return self.surface
        if type_part not in ["sheet", "lower skeletal", "upper skeletal"]:
            raise ValueError(
                f"type_part ({type_part}) must be 'sheet', \
                    'lower skeletal', 'upper skeletal' or 'surface'"
            )
        if self.density is not None:
            self._compute_offset_to_fit_density(
                part_type=type_part,
                resolution=algo_resolution,
            )
        polydata = getattr(self, f"vtk_{type_part.replace(' ', '_')}")()
        polydata = polydata.clean().triangulate()
        density = polydata.volume / self.grid.volume
        logging.info(f"TPMS density = {density:.2%}")

        polydata = rotatePvEuler(
            polydata,
            center=(0, 0, 0),
            psi=self.orientation[0],
            theta=self.orientation[1],
            phi=self.orientation[2],
        )
        return polydata.translate(xyz=self.center)


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
        density: Optional[float] = None,
    ):
        """
        Directions of cell_size and repeat_cell must be taken as the cylindrical \
            coordinate system $\\left(\\rho, \\theta, z\\right)$.

        The $\\theta$ component of cell_size is automatically updated to the \
            closest value that matches the cylindrical periodicity of the TPMS.
        If the $\\theta$ component of repeat_cell is 0 or greater than the \
            periodicity of the TPMS, it is automatically set the correct number \
                to make the full cylinder.

        :param radius: radius of the cylinder on which the center of the TPMS is located
        :param surface_function: tpms function or custom function (f(x, y, z) = 0)
        :param offset: offset of the isosurface to generate thickness
        :param phase_shift: phase shift of the tpms function \
            $f(x + \\phi_x, y + \\phi_y, z + \\phi_z) = 0$
        :param cell_size: float or list of float for each dimension to set unit cell dimensions
        :param repeat_cell: integer or list of integers to repeat the geometry in each dimension
        :param center: center of the geometry
        :param orientation: orientation of the geometry
        :param resolution: unit cell resolution of the grid to compute tpms scalar fields
        """
        self._init_cell_parameters(cell_size, repeat_cell)

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

        super().__init__(
            surface_function=surface_function,
            offset=offset,
            phase_shift=phase_shift,
            cell_size=self.cell_size,
            repeat_cell=self.repeat_cell,
            resolution=resolution,
            density=density,
            center=center,
            orientation=orientation,
        )

    def _create_grid(
        self, x: np.ndarray, y: np.ndarray, z: np.ndarray
    ) -> pv.StructuredGrid:
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
        density: Optional[float] = None,
    ):
        """
        Directions of cell_size and repeat_cell must be taken as the spherical \
            coordinate system $\\left(r, \\theta, \\phi\\right)$.

        The $\\theta$ and $\\phi$ components of cell_size are automatically \
            updated to the closest values that matches the spherical periodicity of the TPMS.
        If the $\\theta$ or $\\phi$ components of repeat_cell are 0 or greater \
            than the periodicity of the TPMS, they are automatically set the correct \
                number to make the full sphere.

        :param radius: radius of the sphere on which the center of the TPMS is located
        :param surface_function: tpms function or custom function (f(x, y, z) = 0)
        :param offset: offset of the isosurface to generate thickness
        :param phase_shift: phase shift of the tpms function \
            $f(x + \\phi_x, y + \\phi_y, z + \\phi_z) = 0$
        :param cell_size: float or list of float for each dimension to set unit cell dimensions
        :param repeat_cell: integer or list of integers to repeat the geometry in each dimension
        :param center: center of the geometry
        :param orientation: orientation of the geometry
        :param resolution: unit cell resolution of the grid to compute tpms scalar fields
        """
        self._init_cell_parameters(cell_size, repeat_cell)

        self.sphere_radius = radius

        unit_theta = self.cell_size[1] / radius
        n_repeat_theta_to_join = int(np.pi / unit_theta)
        self.unit_theta = np.pi / n_repeat_theta_to_join
        self.cell_size[1] = self.unit_theta * radius  # true only on theta = pi/2
        if self.repeat_cell[1] == 0 or self.repeat_cell[1] > n_repeat_theta_to_join:
            logging.info("%d cells repeated in theta direction", n_repeat_theta_to_join)
            self.repeat_cell[1] = n_repeat_theta_to_join

        unit_phi = self.cell_size[2] / radius
        n_repeat_phi_to_join = int(2 * np.pi / unit_phi)
        self.unit_phi = 2 * np.pi / n_repeat_phi_to_join
        self.cell_size[2] = self.unit_phi * radius
        if self.repeat_cell[2] == 0 or self.repeat_cell[2] > n_repeat_phi_to_join:
            logging.info("%d cells repeated in phi direction", n_repeat_phi_to_join)
            self.repeat_cell[2] = n_repeat_phi_to_join

        super().__init__(
            surface_function=surface_function,
            offset=offset,
            phase_shift=phase_shift,
            cell_size=self.cell_size,
            repeat_cell=self.repeat_cell,
            resolution=resolution,
            density=density,
            center=center,
            orientation=orientation,
        )

    def _create_grid(
        self, x: np.ndarray, y: np.ndarray, z: np.ndarray
    ) -> pv.StructuredGrid:
        rho = x + self.sphere_radius
        theta = y * self.unit_theta + np.pi / 2.0
        phi = z * self.unit_phi

        return pv.StructuredGrid(
            rho * np.sin(theta) * np.cos(phi),
            rho * np.sin(theta) * np.sin(phi),
            rho * np.cos(theta),
        )
