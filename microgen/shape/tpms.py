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
from typing import Any, Callable, Literal, Sequence

import cadquery as cq
import numpy as np
import pyvista as pv
from scipy.optimize import root_scalar

from microgen.operations import fuseShapes, rotateEuler, rotatePvEuler

from .basic_geometry import BasicGeometry

logging.basicConfig(level=logging.INFO)
Field = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]

_DIM = 3


class Tpms(BasicGeometry):
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
        offset: float | Field = 0.0,
        phase_shift: Sequence[float] = (0.0, 0.0, 0.0),
        cell_size: float | Sequence[float] = 1.0,
        repeat_cell: int | Sequence[int] = 1,
        resolution: int = 20,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        density: float | None = None,
    ) -> None:
        r"""Class used to generate TPMS geometries (sheet or skeletals parts).

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
            raise DensityError(density)
        self.density = density

    @classmethod
    def offset_from_density(
        cls: type[Tpms],
        surface_function: Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray],
        part_type: Literal["sheet", "lower skeletal", "upper skeletal"],
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
        if self.density is None:
            raise DensityError(self.density)

        if self.density == 1.0:
            offset = (
                self.offset_lim["sheet"][1]
                if part_type == "sheet"
                else self.offset_lim["skeletal"][0]
            )
            self._update_offset(offset)
            return offset

        temp_tpms = Tpms(
            surface_function=self.surface_function,
            resolution=resolution if resolution is not None else self.resolution,
        )
        polydata_func = getattr(temp_tpms, f"vtk_{part_type.replace(' ', '_')}")

        def density(offset: float) -> float:
            temp_tpms._update_offset(offset)  # noqa: SLF001
            return abs(polydata_func().volume)

        part = "skeletal" if "skeletal" in part_type else part_type
        computed_offset = root_scalar(
            lambda offset: density(offset) - self.density,
            bracket=self.offset_lim[part],
        ).root
        self._update_offset(computed_offset)
        return computed_offset

    def _init_cell_parameters(
        self: Tpms,
        cell_size: float | Sequence[float],
        repeat_cell: int | Sequence[int],
    ) -> None:
        if isinstance(cell_size, (float, int)):
            self.cell_size = np.array([cell_size, cell_size, cell_size])
        elif len(cell_size) == _DIM:
            self.cell_size = np.array(cell_size)
        else:
            raise SequenceLengthError(sequence="cell_size", variable_type=float)

        if isinstance(repeat_cell, int):
            self.repeat_cell = np.array([repeat_cell, repeat_cell, repeat_cell])
        elif len(repeat_cell) == _DIM:
            self.repeat_cell = np.array(repeat_cell)
        else:
            raise SequenceLengthError(sequence="repeat_cell", variable_type=int)

    def vtk_sheet(self: Tpms) -> pv.PolyData:
        """Return sheet part."""
        return self.grid.clip_scalar(scalars="lower_surface", invert=False).clip_scalar(
            scalars="upper_surface",
        )

    def vtk_upper_skeletal(self: Tpms) -> pv.PolyData:
        """Return upper skeletal part."""
        return self.grid.clip_scalar(scalars="upper_surface", invert=False)

    def vtk_lower_skeletal(self: Tpms) -> pv.PolyData:
        """Return lower skeletal part."""
        return self.grid.clip_scalar(scalars="lower_surface")

    @property
    def sheet(self: Tpms) -> pv.PolyData:
        """Returns sheet part."""
        if self._sheet is not None:
            return self._sheet

        if self.density is not None:
            self._compute_offset_to_fit_density(part_type="sheet")

        self._sheet = self.vtk_sheet().clean().triangulate()
        return self._sheet

    @property
    def upper_skeletal(self: Tpms) -> pv.PolyData:
        """Returns upper skeletal part."""
        if self._upper_skeletal is not None:
            return self._upper_skeletal

        if self.density is not None:
            self._compute_offset_to_fit_density(part_type="upper skeletal")

        self._upper_skeletal = self.vtk_upper_skeletal().clean().triangulate()
        return self._upper_skeletal

    @property
    def lower_skeletal(self: Tpms) -> pv.PolyData:
        """Returns lower skeletal part."""
        if self._lower_skeletal is not None:
            return self._lower_skeletal

        if self.density is not None:
            self._compute_offset_to_fit_density(part_type="lower skeletal")

        self._lower_skeletal = self.vtk_lower_skeletal().clean().triangulate()
        return self._lower_skeletal

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
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
    ) -> pv.StructuredGrid:
        return pv.StructuredGrid(x, y, z)

    def _compute_tpms_field(self: Tpms) -> None:
        linspaces: list[np.ndarray] = [
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
        tpms_field = self.surface_function(
            k_x * (x + self.phase_shift[0]),
            k_y * (y + self.phase_shift[1]),
            k_z * (z + self.phase_shift[2]),
        )

        self.grid["surface"] = tpms_field.ravel(order="F")
        self._update_offset(self.offset)

    def _update_offset(self: Tpms, offset: float | Callable) -> None:
        if isinstance(offset, float):
            self.offset = offset
        elif isinstance(offset, Callable):
            self.offset = offset(self.grid.x, self.grid.y, self.grid.z).ravel("F")

        self.grid["lower_surface"] = self.grid["surface"] + 0.5 * self.offset
        self.grid["upper_surface"] = self.grid["surface"] - 0.5 * self.offset

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

        return cq.Shell.makeShell(faces)

    def _create_surface(
        self: Tpms,
        isovalue: float | np.ndarray = 0.0,
        smoothing: int = 0,
    ) -> cq.Shell:
        if isinstance(isovalue, (int, float)):
            scalars = self.grid["surface"] - isovalue
        elif isinstance(isovalue, np.ndarray):
            scalars = self.grid["surface"] - isovalue.ravel(order="F")

        mesh = self.grid.contour(isosurfaces=[0.0], scalars=scalars)
        mesh.smooth(n_iter=smoothing, feature_smoothing=True, inplace=True)
        mesh.clean(inplace=True)

        try:
            shell = self._create_shell(mesh=mesh)
        except ValueError as exc:
            raise CreateShellError from exc
        return shell

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
        eps: float,
        smoothing: int,
    ) -> tuple[cq.Shape, cq.Shape]:
        isovalues = [
            -self.offset / 2.0,
            -self.offset / 2.0 + eps,
            self.offset / 2.0 - eps,
            self.offset / 2.0,
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
        eps: float,
        smoothing: int,
    ) -> tuple[cq.Shape, cq.Shape]:
        isovalues = [
            -self.offset / 2.0,
            -self.offset / 2.0 - eps,
        ]

        shells = self._create_surfaces(isovalues, smoothing)

        surface = shells[0]
        test_surface = shells[1]
        return surface, test_surface

    def _generate_upper_skeletal_surfaces(
        self: Tpms,
        eps: float,
        smoothing: int,
    ) -> tuple[cq.Shape, cq.Shape]:
        isovalues = [
            self.offset / 2.0,
            self.offset / 2.0 + eps,
        ]

        shells = self._create_surfaces(isovalues, smoothing)

        surface = shells[0]
        test_surface = shells[1]
        return surface, test_surface

    def _extract_part_from_box(
        self: Tpms,
        type_part: Literal["sheet", "lower skeletal", "upper skeletal", "surface"],
        eps: float,
        smoothing: int,
    ) -> cq.Shape:
        box = cq.Workplane("front").box(*(self.cell_size * self.repeat_cell))

        surface, test_surface = getattr(
            self,
            f"_generate_{type_part.replace(' ', '_')}_surfaces",
        )(eps, smoothing)

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
        type_part: Literal["sheet", "lower skeletal", "upper skeletal", "surface"],
    ) -> None:
        if "skeletal" in type_part:
            if (
                isinstance(self.offset, (int, float)) and self.offset < 0.0
            ):  # scalar offset = 0 is working
                raise NegativeOffsetNotImplementedError
            if isinstance(self.offset, np.ndarray) and np.any(self.offset <= 0.0):
                raise NegativeOffsetNotImplementedError
        elif type_part == "sheet":
            if np.any(self.offset <= 0.0):
                if np.all(self.offset <= 0.0):
                    raise OffsetRangeError(
                        part_type=type_part,
                        offset_bounds=self.offset_lim["sheet"],
                    )
                raise NegativeOffsetNotImplementedError

        part = "skeletal" if "skeletal" in type_part else type_part
        if np.all(self.offset > self.offset_lim[part][1]):
            raise OffsetRangeError(
                part_type=type_part,
                offset_bounds=self.offset_lim[part],
            )

    def generate(
        self: Tpms,
        type_part: Literal[
            "sheet",
            "lower skeletal",
            "upper skeletal",
            "surface",
        ] = "sheet",
        smoothing: int = 0,
        algo_resolution: int | None = None,
        **_: dict[str, Any],
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
            raise TypePartError(type_part)

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

        eps = self.offset / 3.0
        if isinstance(self.offset, float) and self.offset == 0.0:
            eps = 0.1 * np.min(self.cell_size)

        shape = self._extract_part_from_box(type_part, eps, smoothing)

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
        type_part: Literal[
            "sheet",
            "lower skeletal",
            "upper skeletal",
            "surface",
        ] = "sheet",
        algo_resolution: int | None = None,
        **_: dict[str, Any],
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
            raise TypePartError(type_part)
        if self.density is not None:
            self._compute_offset_to_fit_density(
                part_type=type_part,
                resolution=algo_resolution,
            )
        polydata = getattr(self, f"vtk_{type_part.replace(' ', '_')}")()
        polydata = polydata.clean().triangulate()

        polydata = rotatePvEuler(
            polydata,
            center=(0, 0, 0),
            psi=self.orientation[0],
            theta=self.orientation[1],
            phi=self.orientation[2],
        )
        return polydata.translate(xyz=self.center)

    def generateVtk(  # noqa: N802
        self: Tpms,
        type_part: Literal[
            "sheet",
            "lower skeletal",
            "upper skeletal",
            "surface",
        ] = "sheet",
        algo_resolution: int | None = None,
        **_: dict[str, Any],
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
        offset: float | Field = 0.0,
        phase_shift: Sequence[float] = (0.0, 0.0, 0.0),
        cell_size: float | Sequence[float] = 1.0,
        repeat_cell: int | Sequence[int] = 1,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
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
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
    ) -> pv.StructuredGrid:
        rho = x + self.cylinder_radius
        theta = y * self.unit_theta

        return pv.StructuredGrid(rho * np.cos(theta), rho * np.sin(theta), z)


class SphericalTpms(Tpms):
    """Class used to generate spherical TPMS geometries (sheet or skeletals parts)."""

    def __init__(  # noqa: PLR0913
        self: SphericalTpms,
        radius: float,
        surface_function: Field,
        offset: float | Field = 0.0,
        phase_shift: Sequence[float] = (0.0, 0.0, 0.0),
        cell_size: float | Sequence[float] = 1.0,
        repeat_cell: int | Sequence[int] = 1,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
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
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
    ) -> pv.StructuredGrid:
        rho = x + self.sphere_radius
        theta = y * self.unit_theta + np.pi / 2.0
        phi = z * self.unit_phi

        return pv.StructuredGrid(
            rho * np.sin(theta) * np.cos(phi),
            rho * np.sin(theta) * np.sin(phi),
            rho * np.cos(theta),
        )


class DensityError(ValueError):
    """Raised when the density is not between 0 and 1."""

    def __init__(self: DensityError, density: float | None) -> None:
        """Initialize the error message."""
        self.message = f"density must be between 0 and 1. Given: {density}"
        super().__init__(self.message)


class SequenceLengthError(ValueError):
    """Raised when the length of the sequence is not 3."""

    def __init__(
        self: SequenceLengthError,
        sequence: Sequence,
        variable_type: type,
    ) -> None:
        """Initialize the error message."""
        self.message = (
            f"Sequence must have a length of 3 {variable_type}. Given: {sequence}"
        )
        super().__init__(self.message)


class CreateShellError(ValueError):
    """Raised when the shell cannot be created."""

    def __init__(self: CreateShellError) -> None:
        """Initialize the error message."""
        self.message = "Cannot create shell, try to use a higher smoothing value."
        super().__init__(self.message)


class NegativeOffsetNotImplementedError(NotImplementedError):
    """Raised when the offset is negative."""

    def __init__(self: NegativeOffsetNotImplementedError) -> None:
        """Initialize the error message."""
        self.message = "generating part with a negative or zero offset\
              value is not implemented yet"
        super().__init__(self.message)


class TypePartError(ValueError):
    """Raised when the type_part is not valid."""

    def __init__(self: TypePartError, type_part: str) -> None:
        """Initialize the error message."""
        self.message = f"type_part ({type_part}) must be 'sheet', 'lower skeletal',\
              'upper skeletal' or 'surface'"
        super().__init__(self.message)


class OffsetRangeError(ValueError):
    """Raised when the offset is not valid for the sheet part."""

    def __init__(
        self: OffsetRangeError,
        part_type: str,
        offset_bounds: tuple[float, float],
    ) -> None:
        """Initialize the error message."""
        self.message = f"offset must be greater than {offset_bounds[0]} to \
            generate '{part_type}' part and lower than {offset_bounds[1]}"
        super().__init__(self.message)
