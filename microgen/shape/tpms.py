"""TPMS.

=============================================
TPMS (:mod:`microgen.shape.tpms`)
=============================================

.. jupyter-execute::
   :hide-code:

   pyvista.global_theme.smooth_shading = True
   pyvista.global_theme.split_sharp_edges = True

"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Callable, Literal, Sequence

import cadquery as cq
import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.optimize import root_scalar

from microgen.operations import fuseShapes, rotateEuler, rotatePvEuler

from .shape import Shape

if TYPE_CHECKING:
    from microgen.shape import KwargsGenerateType, TpmsPartType, Vector3DType
from .tpms_grading import OffsetGrading

logging.basicConfig(level=logging.INFO)
Field = Callable[
    [npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]],
    npt.NDArray[np.float64],
]

_DIM = 3


class Tpms(Shape):
    """Triply Periodical Minimal Surfaces.

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

    def __init__(  # noqa: PLR0913
        self: Tpms,
        surface_function: Field,
        offset: float | OffsetGrading | Field | None = None,
        phase_shift: Sequence[float] = (0.0, 0.0, 0.0),
        cell_size: float | Sequence[float] = 1.0,
        repeat_cell: int | Sequence[int] = 1,
        resolution: int = 20,
        density: float | None = None,
        **kwargs: Vector3DType,
    ) -> None:
        r"""Class used to generate TPMS geometries (sheet or skeletals parts).

        TPMS are created by default in a cube.
        The geometry of the cube can be modified using 'cell_size' parameter.
        The number of repetitions in each direction of the created geometry \
            can be modified with the 'repeat_cell' parameter.

        :param surface_function: tpms function or custom function (f(x, y, z) = 0)
        :param offset: offset of the isosurface to generate thickness
        :param phase_shift: phase shift of the isosurface \
            $f(x + \\phi_x, y + \\phi_y, z + \\phi_z, t) = 0$
        :param cell_size: float or list of float for each dimension to set unit\
              cell dimensions
        :param repeat_cell: integer or list of integers to repeat the geometry\
              in each dimension
        :param resolution: unit cell resolution of the grid to compute tpms\
              scalar fields
        :param density: density percentage of the generated geometry (0 < density < 1) \
            If density is given, the offset is automatically computed to fit the\
                  density (performance is slower than when using the offset)
        """
        if offset is not None and density is not None:
            err_msg = (
                "offset and density cannot be given at the same time. Give only one."
            )
            raise ValueError(err_msg)
        if offset is None and density is None:
            err_msg = "offset or density must be given. Give one of them."
            raise ValueError(err_msg)

        super().__init__(**kwargs)

        self.surface_function = surface_function
        self._offset = offset if offset is not None else 0.0
        self.phase_shift = phase_shift

        self.grid: pv.StructuredGrid
        self._grid_sheet: pv.UnstructuredGrid = None
        self._grid_upper_skeletal: pv.UnstructuredGrid = None
        self._grid_lower_skeletal: pv.UnstructuredGrid = None
        self._surface: pv.PolyData = None

        self._init_cell_parameters(cell_size, repeat_cell)

        self.resolution = resolution

        self._compute_tpms_field()

        min_field = np.min(self.grid["surface"])
        max_field = np.max(self.grid["surface"])
        self.offset_lim = {
            "sheet": (
                0.0,
                2.0 * max(-min_field, max_field),
            ),
            "skeletal": (
                2.0 * min_field,
                2.0 * max_field,
            ),
        }

        if density is not None and not 0.0 < density <= 1.0:
            err_msg = f"density must be between 0 and 1. Given: {density}"
            raise ValueError(err_msg)
        self.density = density

        self._offset: float | npt.NDArray | None
        if density is not None:
            self._offset = None
        else:
            self.offset = offset  # call setter
        self.offset_updated: bool

    @classmethod
    def offset_from_density(
        cls: type[Tpms],
        surface_function: Field,
        part_type: TpmsPartType,
        density: float,
        resolution: int = 20,
    ) -> float:
        """Return the offset corresponding to the required density.

        :param surface_function: tpms function
        :param part_type: type of the part (sheet, lower skeletal or upper skeletal)
        :param density: Required density, 0.5 for 50%
        :param resolution: resolution of the tpms used to compute the offset

        :return: corresponding offset value
        """
        return Tpms(  # noqa: SLF001
            surface_function=surface_function,
            density=density,
        )._compute_offset_to_fit_density(part_type=part_type, resolution=resolution)

    def _compute_offset_to_fit_density(
        self: Tpms,
        part_type: Literal["sheet", "lower skeletal", "upper skeletal"],
        resolution: int | None = None,
    ) -> float:
        """Compute the offset to fit the required density."""
        if self.density is None:
            err_msg = f"density must be between 0 and 1. Given: {self.density}"
            raise ValueError(err_msg)

        if self.density == 1.0:
            self.offset = (
                self.offset_lim["sheet"][1]
                if part_type == "sheet"
                else self.offset_lim["skeletal"][0]
            )
            return self.offset

        temp_tpms = Tpms(
            surface_function=self.surface_function,
            offset=0.0,
            resolution=resolution if resolution is not None else self.resolution,
        )

        def density(offset: float) -> float:
            temp_tpms.offset = offset
            grid_part = getattr(temp_tpms, f"grid_{part_type.replace(' ', '_')}")
            return abs(grid_part.volume)

        part = "skeletal" if "skeletal" in part_type else part_type
        computed_offset = root_scalar(
            lambda offset: density(offset) - self.density,
            bracket=self.offset_lim[part],
        ).root
        self.offset = computed_offset
        return computed_offset

    def _init_cell_parameters(
        self: Tpms,
        cell_size: float | Sequence[float],
        repeat_cell: int | Sequence[int],
    ) -> None:
        """Initialize cell_size and repeat_cell so that they are arrays of length 3."""
        if isinstance(cell_size, (float, int)):
            self.cell_size = np.array([cell_size, cell_size, cell_size])
        elif len(cell_size) == _DIM:
            self.cell_size = np.array(cell_size)
        else:
            err_msg = f"`cell_size` must have a length of 3 floats. Given: {cell_size}"
            raise ValueError(err_msg)

        if isinstance(repeat_cell, int):
            self.repeat_cell = np.array([repeat_cell, repeat_cell, repeat_cell])
        elif len(repeat_cell) == _DIM:
            self.repeat_cell = np.array(repeat_cell)
        else:
            err_msg = (
                "`repeat_cell` must have a length of 3 integers. "
                f"Given: {repeat_cell}"
            )
            raise ValueError(err_msg)

    @property
    def grid_sheet(self: Tpms) -> pv.UnstructuredGrid:
        """Return sheet part."""
        if self._grid_sheet is not None and not self.offset_updated:
            return self._grid_sheet
        self.offset_updated = False

        if self.density is not None:
            self._compute_offset_to_fit_density(part_type="sheet")

        self._grid_sheet = self.grid.clip_scalar(
            scalars="lower_surface",
            invert=False,
        ).clip_scalar(
            scalars="upper_surface",
        )
        return self._grid_sheet

    @property
    def grid_upper_skeletal(self: Tpms) -> pv.UnstructuredGrid:
        """Return upper skeletal part."""
        if self._grid_upper_skeletal is not None and not self.offset_updated:
            return self._grid_upper_skeletal
        self.offset_updated = False

        if self.density is not None:
            self._compute_offset_to_fit_density(part_type="upper skeletal")

        self._grid_upper_skeletal = self.grid.clip_scalar(
            scalars="upper_surface",
            invert=False,
        )
        return self._grid_upper_skeletal

    @property
    def grid_lower_skeletal(self: Tpms) -> pv.UnstructuredGrid:
        """Return lower skeletal part."""
        if self._grid_lower_skeletal is not None and not self.offset_updated:
            return self._grid_lower_skeletal
        self.offset_updated = False

        if self.density is not None:
            self._compute_offset_to_fit_density(part_type="lower skeletal")

        self._grid_lower_skeletal = self.grid.clip_scalar(scalars="lower_surface")
        return self._grid_lower_skeletal

    @property
    def sheet(self: Tpms) -> pv.PolyData:
        """Return sheet part."""
        return self.grid_sheet.extract_surface().clean().triangulate()

    @property
    def upper_skeletal(self: Tpms) -> pv.PolyData:
        """Return upper skeletal part."""
        return self.grid_upper_skeletal.extract_surface().clean().triangulate()

    @property
    def lower_skeletal(self: Tpms) -> pv.PolyData:
        """Return lower skeletal part."""
        return self.grid_lower_skeletal.extract_surface().clean().triangulate()

    @property
    def skeletals(self: Tpms) -> tuple[pv.PolyData, pv.PolyData]:
        """Returns both skeletal parts."""
        return (self.upper_skeletal, self.lower_skeletal)

    @property
    def surface(self: Tpms) -> pv.PolyData:
        """Returns isosurface f(x, y, z) = 0."""
        if self._surface is not None:
            return self._surface

        self._surface = self.grid.contour(
            isosurfaces=[0.0],
            scalars="surface",
        ).triangulate()
        return self._surface

    def _create_grid(
        self: Tpms,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> pv.StructuredGrid:
        """Return the structured cartesian grid of the TPMS."""
        grid = pv.StructuredGrid(x, y, z)
        grid["coords"] = np.c_[
            x.ravel(order="F"),
            y.ravel(order="F"),
            z.ravel(order="F"),
        ]
        return grid

    def _compute_tpms_field(self: Tpms) -> None:
        """Compute the TPMS scalar field on the grid."""
        linspaces: list[npt.NDArray[np.float64]] = [
            np.linspace(
                -0.5 * cell_size_axis * repeat_cell_axis,
                0.5 * cell_size_axis * repeat_cell_axis,
                self.resolution * repeat_cell_axis,
            )
            for repeat_cell_axis, cell_size_axis in zip(
                self.repeat_cell,
                self.cell_size,
            )
        ]

        x, y, z = np.meshgrid(*linspaces)

        self.grid = self._create_grid(x, y, z)

        k_x, k_y, k_z = 2.0 * np.pi / self.cell_size
        x, y, z = self.grid["coords"].T
        tpms_field = self.surface_function(
            k_x * (x + self.phase_shift[0]),
            k_y * (y + self.phase_shift[1]),
            k_z * (z + self.phase_shift[2]),
        )

        self.grid["surface"] = tpms_field.ravel(order="F")

    def _update_grid_offset(self: Tpms) -> None:
        self.grid["lower_surface"] = self.grid["surface"] + 0.5 * self.offset
        self.grid["upper_surface"] = self.grid["surface"] - 0.5 * self.offset

    @property
    def offset(self: Tpms) -> float | npt.NDArray[np.float64] | None:
        """Returns the offset value."""
        return self._offset

    @offset.setter
    def offset(
        self: Tpms,
        offset: float | npt.NDArray[np.float64] | OffsetGrading | Field,
    ) -> None:
        if isinstance(offset, (int, float, np.ndarray)):
            self._offset = offset
        elif isinstance(offset, OffsetGrading):
            self._offset = offset.compute_offset(self.grid)
        elif callable(offset):
            self._offset = offset(self.grid.x, self.grid.y, self.grid.z).ravel("F")
        else:
            err_msg = "offset must be a float, a numpy array or a callable"
            raise TypeError(err_msg)

        self._update_grid_offset()
        self.offset_updated = True

    def _create_shell(self: Tpms, mesh: pv.PolyData) -> cq.Shell:
        if not mesh.is_all_triangles:
            mesh.triangulate(inplace=True)  # useless ?
        triangles = mesh.faces.reshape(-1, 4)[:, 1:]
        triangles = np.c_[triangles, triangles[:, 0]]

        faces = []
        for tri in triangles:
            lines = [
                cq.Edge.makeLine(
                    cq.Vector(*mesh.points[start]),
                    cq.Vector(*mesh.points[end]),
                )
                for start, end in zip(tri[:], tri[1:])
            ]

            wire = cq.Wire.assembleEdges(lines)
            faces.append(cq.Face.makeFromWires(wire))

        try:
            shell = cq.Shell.makeShell(faces)
        except ValueError as err:
            err_msg = "Failed to create the shell, \
                try to increase the resolution or the smoothing."
            raise ShellCreationError(err_msg) from err
        return shell

    def _create_surface(
        self: Tpms,
        isovalue: float | npt.NDArray[np.float64] = 0.0,
        smoothing: int = 0,
    ) -> cq.Shell:
        """Create a TPMS surface for the given isovalue."""
        if isinstance(isovalue, (int, float)):
            scalars = self.grid["surface"] - isovalue
        elif isinstance(isovalue, np.ndarray):
            scalars = self.grid["surface"] - isovalue.ravel(order="F")

        mesh = self.grid.contour(isosurfaces=[0.0], scalars=scalars)
        mesh.smooth(n_iter=smoothing, feature_smoothing=True, inplace=True)
        mesh.clean(inplace=True)

        return self._create_shell(mesh=mesh)

    def _create_surfaces(
        self: Tpms,
        isovalues: list[float],
        smoothing: int = 0,
    ) -> list[cq.Shell]:
        """Create TPMS surfaces for the corresponding isovalue.

        :param isovalues: list of isovalues corresponding to the required surfaces
        :param smoothing: smoothing loop iterations
        :param verbose: display progressbar of the conversion to CadQuery object

        :return: list of CadQuery Shell objects of the required TPMS surfaces
        """
        shells = []
        for i, isovalue in enumerate(isovalues):
            logging.info("\nGenerating surface (%d/%d)", i + 1, len(isovalues))
            shell = self._create_surface(
                isovalue=isovalue,
                smoothing=smoothing,
            )
            shells.append(shell)

        return shells

    def _generate_sheet_surfaces(
        self: Tpms,
        smoothing: int,
    ) -> tuple[cq.Shape, cq.Shape]:
        """Generate the surfaces to create the sheet part of the TPMS."""
        isovalues = [
            -0.5 * self.offset,
            -0.25 * self.offset,
            0.25 * self.offset,
            0.5 * self.offset,
        ]

        shells = self._create_surfaces(isovalues, smoothing)

        lower_surface = shells[0]
        lower_test_surface = shells[1]
        upper_test_surface = shells[2]
        upper_surface = shells[3]

        surface = lower_surface.fuse(upper_surface)
        test_surface = lower_test_surface.fuse(upper_test_surface)
        return surface, test_surface

    def _generate_lower_skeletal_surfaces(
        self: Tpms,
        smoothing: int,
    ) -> tuple[cq.Shape, cq.Shape]:
        """Generate the surfaces to create the lower skeletal part of the TPMS."""
        min_offset = 2.0 * np.min(self.grid["surface"])
        isovalues = [
            -0.5 * self.offset,
            -0.25 * (self.offset - min_offset),
        ]

        shells = self._create_surfaces(isovalues, smoothing)

        surface = shells[0]
        test_surface = shells[1]
        return surface, test_surface

    def _generate_upper_skeletal_surfaces(
        self: Tpms,
        smoothing: int,
    ) -> tuple[cq.Shape, cq.Shape]:
        """Generate the surfaces to create the upper skeletal part of the TPMS."""
        min_offset = 2.0 * np.min(self.grid["surface"])
        isovalues = [
            0.5 * self.offset,
            0.25 * (self.offset - min_offset),
        ]

        shells = self._create_surfaces(isovalues, smoothing)

        surface = shells[0]
        test_surface = shells[1]
        return surface, test_surface

    def _extract_part_from_box(
        self: Tpms,
        type_part: TpmsPartType,
        smoothing: int,
    ) -> cq.Shape:
        """Extract the required part from the box."""
        box = cq.Workplane("front").box(*(self.cell_size * self.repeat_cell))

        surface, test_surface = getattr(
            self,
            f"_generate_{type_part.replace(' ', '_')}_surfaces",
        )(smoothing)

        splitted_box = box.split(surface)
        tpms_solids = splitted_box.solids().all()

        # split each solid with the test surface to identify
        # to what part type the solid belongs to
        list_solids = [
            (solid.split(test_surface).solids().size(), solid.val())
            for solid in tpms_solids
        ]

        # if the number of shapes is greater than 1, it means that the solid is split
        # so it belongs to the required part
        part_solids = [solid for (number, solid) in list_solids if number > 1]
        part_shapes = [cq.Shape(solid.wrapped) for solid in part_solids]
        return fuseShapes(
            cqShapeList=part_shapes,
            retain_edges=False,  # True or False ?
        )

    def _check_offset(
        self: Tpms,
        type_part: TpmsPartType,
    ) -> None:
        """Check if the offset is valid for the required part."""
        if "skeletal" in type_part:
            if (
                isinstance(self.offset, (int, float)) and self.offset < 0.0
            ):  # scalar offset = 0 is working
                err_msg = "generating part with a negative or zero offset\
                    value is not implemented yet"
                raise NotImplementedError(err_msg)
            if isinstance(self.offset, np.ndarray) and np.any(self.offset <= 0.0):
                err_msg = "generating part with a negative or zero offset\
                    value is not implemented yet"
                raise NotImplementedError(err_msg)
        elif type_part == "sheet":
            if np.any(self.offset <= 0.0):
                if np.all(self.offset <= 0.0):
                    err_msg = (
                        f"offset must be greater than "
                        f"{self.offset_lim[type_part][0]} to generate '{type_part}' "
                        f"part and lower than {self.offset_lim[type_part][1]}"
                    )
                    raise ValueError(err_msg)
                err_msg = "generating part with a negative or zero offset\
                    value is not implemented yet"
                raise NotImplementedError(err_msg)

        part = "skeletal" if "skeletal" in type_part else type_part
        if np.all(self.offset > self.offset_lim[part][1]):
            err_msg = f"offset must be greater than {self.offset_lim[part][0]} to \
                generate '{type_part}' part and lower than {self.offset_lim[part][1]}"
            raise ValueError(err_msg)

    def generate(
        self: Tpms,
        type_part: TpmsPartType = "sheet",
        smoothing: int = 0,
        algo_resolution: int | None = None,
        **_: KwargsGenerateType,
    ) -> cq.Shape:
        """Generate CadQuery Shape object of the required TPMS part.

        :param type_part: part of the TPMS desired \
            ('sheet', 'lower skeletal', 'upper skeletal' or 'surface')
        :param smoothing: smoothing loop iterations
        :param verbose: display progressbar of the conversion to CadQuery object
        :param algo_resolution: if offset must be computed to fit density, \
            resolution of the temporary TPMS used to compute the offset

        :return: CadQuery Shape object of the required TPMS part
        """
        if type_part not in ["sheet", "lower skeletal", "upper skeletal", "surface"]:
            err_msg = (
                f"type_part ({type_part}) must be 'sheet', 'lower skeletal', "
                "'upper skeletal' or 'surface'"
            )
            raise ValueError(err_msg)

        if type_part == "surface":
            if self.offset != 0.0:
                logging.warning("offset is ignored for 'surface' part")
            if self.density is not None:
                logging.warning("density is ignored for 'surface' part")
            return self._create_surface(
                isovalue=0,
                smoothing=smoothing,
            )

        if self.density is not None:
            self._compute_offset_to_fit_density(
                part_type=type_part,
                resolution=algo_resolution,
            )

        self._check_offset(type_part)

        shape = self._extract_part_from_box(type_part, smoothing)

        shape = rotateEuler(
            obj=shape,
            center=(0, 0, 0),
            psi=self.orientation[0],
            theta=self.orientation[1],
            phi=self.orientation[2],
        )
        return shape.translate(self.center)

    def generate_vtk(
        self: Tpms,
        type_part: TpmsPartType = "sheet",
        algo_resolution: int | None = None,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate VTK PolyData object of the required TPMS part.

        :param type_part: part of the TPMS desireds
        :param algo_resolution: if offset must be computed to fit density, \
            resolution of the temporary TPMS used to compute the offset

        :return: VTK PolyData object of the required TPMS part
        """
        if type_part == "surface":
            return self.surface
        if type_part not in ["sheet", "lower skeletal", "upper skeletal"]:
            err_msg = (
                f"type_part ({type_part}) must be 'sheet', 'lower skeletal', "
                "'upper skeletal' or 'surface'"
            )
            raise ValueError(err_msg)
        if self.density is not None:
            self._compute_offset_to_fit_density(
                part_type=type_part,
                resolution=algo_resolution,
            )
        polydata = getattr(self, type_part.replace(" ", "_"))

        polydata = rotatePvEuler(
            polydata,
            center=(0, 0, 0),
            psi=self.orientation[0],
            theta=self.orientation[1],
            phi=self.orientation[2],
        )
        return polydata.translate(xyz=self.center)

    def generate_grid_vtk(
        self: Tpms,
        type_part: TpmsPartType = "sheet",
        algo_resolution: int | None = None,
        **_: KwargsGenerateType,
    ) -> pv.UnstructuredGrid:
        """Generate VTK UnstructuredGrid object of the required TPMS part."""
        if type_part not in ["sheet", "lower skeletal", "upper skeletal"]:
            err_msg = (
                f"type_part ({type_part}) must be 'sheet', 'lower skeletal', "
                "'upper skeletal'"
            )
            raise ValueError(err_msg)

        if self.density is not None:
            self._compute_offset_to_fit_density(
                part_type=type_part,
                resolution=algo_resolution,
            )

        return getattr(self, f"grid_{type_part.replace(' ', '_')}")

    def generateVtk(  # noqa: N802
        self: Tpms,
        type_part: TpmsPartType = "sheet",
        algo_resolution: int | None = None,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk(
            type_part=type_part,
            algo_resolution=algo_resolution,
        )


class CylindricalTpms(Tpms):
    """Class used to generate cylindrical TPMS geometries (sheet or skeletals parts)."""

    def __init__(  # noqa: PLR0913
        self: CylindricalTpms,
        radius: float,
        surface_function: Field,
        offset: float | OffsetGrading | Field | None = None,
        phase_shift: Sequence[float] = (0.0, 0.0, 0.0),
        cell_size: float | Sequence[float] = 1.0,
        repeat_cell: int | Sequence[int] = 1,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType = (0, 0, 0),
        resolution: int = 20,
        density: float | None = None,
    ) -> None:
        r"""Cylindrical TPMS geometry.

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
        :param cell_size: float or list of float for each dimension to\
              set unit cell dimensions
        :param repeat_cell: integer or list of integers to repeat the\
              geometry in each dimension
        :param center: center of the geometry
        :param orientation: orientation of the geometry
        :param resolution: unit cell resolution of the grid to compute\
              tpms scalar fields
        """
        self._init_cell_parameters(cell_size, repeat_cell)

        self.cylinder_radius = radius

        unit_theta = self.cell_size[1] / radius
        n_repeat_to_full_circle = int(2 * np.pi / unit_theta)
        self.unit_theta = 2 * np.pi / n_repeat_to_full_circle
        self.cell_size[1] = self.unit_theta * radius
        if self.repeat_cell[1] == 0 or self.repeat_cell[1] > n_repeat_to_full_circle:
            logging.info(
                "%d cells repeated in circular direction",
                n_repeat_to_full_circle,
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
        self: CylindricalTpms,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> pv.StructuredGrid:
        """Return the structured cylindrical grid of the TPMS."""
        rho = x + self.cylinder_radius
        theta = y * self.unit_theta

        grid = pv.StructuredGrid(rho * np.cos(theta), rho * np.sin(theta), z)
        grid["coords"] = np.c_[
            x.ravel(order="F"),
            y.ravel(order="F"),
            z.ravel(order="F"),
        ]
        return grid


class SphericalTpms(Tpms):
    """Class used to generate spherical TPMS geometries (sheet or skeletals parts)."""

    def __init__(  # noqa: PLR0913
        self: SphericalTpms,
        radius: float,
        surface_function: Field,
        offset: float | OffsetGrading | Field | None = None,
        phase_shift: Sequence[float] = (0.0, 0.0, 0.0),
        cell_size: float | Sequence[float] = 1.0,
        repeat_cell: int | Sequence[int] = 1,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType = (0, 0, 0),
        resolution: int = 20,
        density: float | None = None,
    ) -> None:
        r"""Spherical TPMS geometry.

        Directions of cell_size and repeat_cell must be taken as the spherical \
            coordinate system $\\left(r, \\theta, \\phi\\right)$.

        The $\\theta$ and $\\phi$ components of cell_size are automatically \
            updated to the closest values that matches the spherical periodicity\
                  of the TPMS.
        If the $\\theta$ or $\\phi$ components of repeat_cell are 0 or greater \
            than the periodicity of the TPMS, they are automatically set the correct \
                number to make the full sphere.

        :param radius: radius of the sphere on which the center of the TPMS is located
        :param surface_function: tpms function or custom function (f(x, y, z) = 0)
        :param offset: offset of the isosurface to generate thickness
        :param phase_shift: phase shift of the tpms function \
            $f(x + \\phi_x, y + \\phi_y, z + \\phi_z) = 0$
        :param cell_size: float or list of float for each dimension\
              to set unit cell dimensions
        :param repeat_cell: integer or list of integers to repeat the\
              geometry in each dimension
        :param center: center of the geometry
        :param orientation: orientation of the geometry
        :param resolution: unit cell resolution of the grid to compute\
              tpms scalar fields
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
        self: SphericalTpms,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> pv.StructuredGrid:
        """Return the structured spherical grid of the TPMS."""
        rho = x + self.sphere_radius
        theta = y * self.unit_theta + np.pi / 2.0
        phi = z * self.unit_phi

        grid = pv.StructuredGrid(
            rho * np.sin(theta) * np.cos(phi),
            rho * np.sin(theta) * np.sin(phi),
            rho * np.cos(theta),
        )
        grid["coords"] = np.c_[
            x.ravel(order="F"),
            y.ravel(order="F"),
            z.ravel(order="F"),
        ]
        return grid


class Infill(Tpms):
    """Generate a TPMS infill inside a given object."""

    def __init__(  # noqa: PLR0913
        self: Infill,
        obj: pv.PolyData,
        surface_function: Field,
        offset: float | OffsetGrading | Field | None = None,
        cell_size: float | Sequence[float] | npt.NDArray[np.float64] | None = None,
        repeat_cell: int | Sequence[int] | npt.NDArray[np.int8] | None = None,
        phase_shift: Sequence[float] = (0.0, 0.0, 0.0),
        resolution: int = 20,
        density: float | None = None,
    ) -> None:
        r"""Initialize the Infill object.

        :param obj: object in which the infill is generated. Normals must be oriented\
                towards the outside of the object. Use the `flip_normals` method if\
                     needed.
        :param surface_function: tpms function or custom function (f(x, y, z) = 0)
        :param offset: offset of the isosurface to generate thickness
        :param cell_size: float or list of float for each dimension to set\
              unit cell dimensions
        :param repeat_cell: integer or list of integers to repeat the geometry\
              in each dimension
        :param phase_shift: phase shift of the tpms function \
            $f(x + \\phi_x, y + \\phi_y, z + \\phi_z) = 0$
        :param resolution: unit cell resolution of the grid to compute tpms scalar\
              fields
        :param density: density percentage of the generated geometry (0 < density < 1) \
            If density is given, the offset is automatically computed to fit the\
                  density (performance is slower than when using the offset)
        """
        self.obj = obj
        bounds = np.array(obj.bounds)

        margin_factor = 1.001  # to avoid the object surface that can create issues
        obj_dim = margin_factor * (bounds[1::2] - bounds[::2])  # [dim_x, dim_y, dim_z]

        if cell_size is not None and repeat_cell is not None:
            err_msg = (
                "cell_size and repeat_cell cannot be given at the same time, "
                "one is computed from the other."
            )
            raise ValueError(err_msg)

        if cell_size is not None:
            repeat_cell = np.round(obj_dim / cell_size).astype(int)
        elif repeat_cell is not None:
            cell_size = obj_dim / repeat_cell

        if np.any(cell_size > obj_dim):
            err_msg = (
                "cell_size must be lower than the object dimensions. "
                f"Given: {cell_size}, Object dimensions: {obj_dim}"
            )
            raise ValueError(err_msg)

        self._init_cell_parameters(cell_size, repeat_cell)
        super().__init__(
            surface_function=surface_function,
            offset=offset,
            phase_shift=phase_shift,
            cell_size=self.cell_size,
            repeat_cell=self.repeat_cell,
            resolution=resolution,
            density=density,
        )

    def _create_grid(
        self: Infill,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> pv.StructuredGrid:
        grid = super()._create_grid(
            x + self.obj.center[0],
            y + self.obj.center[1],
            z + self.obj.center[2],
        )
        grid = grid.clip_surface(self.obj)
        logging.info("Grid resolution: %s points", grid.n_points)
        return grid


class ShellCreationError(Exception):
    """Error raised when the shell creation fails."""

    def __init__(self: ShellCreationError, message: str) -> None:
        """Initialize the ShellCreationError."""
        super().__init__(message)
