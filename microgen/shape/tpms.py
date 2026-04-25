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

from microgen.operations import fuseShapes, rotate

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
        # Stores the offset callable when one is provided, so the F-rep path
        # can re-evaluate variable thickness on its own marching-cubes grid.
        self._offset_func: Field | None = (
            offset if (offset is not None and callable(offset) and not isinstance(offset, OffsetGrading)) else None
        )
        self.phase_shift = phase_shift

        self.grid: pv.StructuredGrid
        self._grid_sheet: pv.UnstructuredGrid = None
        self._grid_upper_skeletal: pv.UnstructuredGrid = None
        self._grid_lower_skeletal: pv.UnstructuredGrid = None
        self._surface: pv.PolyData = None

        self._init_cell_parameters(cell_size, repeat_cell)

        self.resolution = resolution

        self._compute_tpms_field()
        self._setup_frep_field()

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

    def _density_envelope_volume(self: Tpms) -> float:
        """Return the *envelope* volume that ``density`` is measured against.

        For a plain :class:`Tpms` this is the cartesian cell volume.  Subclasses
        with a non-cartesian envelope (:class:`Infill`, :class:`Conformal`)
        override this to use ``obj.volume`` / ``envelope.volume``, so that
        ``density`` consistently means *fraction of the envelope filled by the
        TPMS part* — not fraction of the cartesian bounding box.
        """
        return abs(self.grid.volume)

    def _compute_offset_to_fit_density(
        self: Tpms,
        part_type: Literal["sheet", "lower skeletal", "upper skeletal"],
        resolution: int | None = None,
    ) -> float:
        """Compute the offset that yields the requested density.

        Searches with the same F-rep ``generate_vtk`` pipeline the user
        invokes, so the offset returned actually reproduces the requested
        density at the user-facing resolution.  When the target density is
        too high to reach at this resolution (marching cubes saturates
        slightly below 1.0 due to surface-grid discretization), falls back
        to the offset limit instead of failing the bracket.

        Density is measured against :meth:`_density_envelope_volume`, so that
        :class:`Infill` / :class:`Conformal` density is relative to the input
        envelope volume rather than the cartesian bbox.
        """
        if self.density is None:
            err_msg = f"density must be between 0 and 1. Given: {self.density}"
            raise ValueError(err_msg)

        part = "skeletal" if "skeletal" in part_type else part_type
        if self.density == 1.0:
            self.offset = (
                self.offset_lim["sheet"][1]
                if part_type == "sheet"
                else self.offset_lim["skeletal"][0]
            )
            return self.offset

        # Density is measured directly on ``self`` so subclass envelopes are
        # respected without having to clone the instance.  ``self.density``
        # is cleared during the search so that ``generate_vtk`` does NOT
        # recurse back into this method.  The final ``computed_offset`` is
        # set as the permanent offset before returning (and ``self.density``
        # is restored so that downstream callers can re-query it).
        envelope_volume = self._density_envelope_volume()
        target_density = self.density
        self.density = None

        def density(offset: float) -> float:
            self.offset = offset  # setter; also clears _offset_func
            mesh = self.generate_vtk(type_part=part_type)
            return abs(mesh.volume) / envelope_volume

        bracket = self.offset_lim[part]
        try:
            computed_offset = root_scalar(
                lambda offset: density(offset) - target_density,
                bracket=bracket,
            ).root
        except ValueError:
            # Bracket sign mismatch — usually because target density exceeds
            # the achievable max at this resolution (or the envelope shape).
            # Pick the bracket endpoint that lies on the side of target_density.
            d_lo = density(bracket[0]) - target_density
            d_hi = density(bracket[1]) - target_density
            computed_offset = bracket[1] if abs(d_hi) <= abs(d_lo) else bracket[0]
        finally:
            self.density = target_density

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
                f"`repeat_cell` must have a length of 3 integers. Given: {repeat_cell}"
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

    def _finalize_frep(
        self: Tpms,
        raw_field: Field,
        bounds: tuple[float, float, float, float, float, float],
    ) -> None:
        """Normalize a raw field to SDF and set ``_func`` / ``_bounds``."""
        from .implicit_ops import from_field, normalize_to_sdf  # noqa: PLC0415

        self._raw_field_func = raw_field
        sdf_shape = normalize_to_sdf(from_field(raw_field))
        self._func = sdf_shape.func
        self._bounds = bounds

    def _setup_frep_field(self: Tpms) -> None:
        """Build the F-rep implicit field (SDF-normalized) for this TPMS."""
        k_x, k_y, k_z = 2.0 * np.pi / self.cell_size
        ps = self.phase_shift

        def _raw_field(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            return self.surface_function(
                k_x * (x + ps[0]),
                k_y * (y + ps[1]),
                k_z * (z + ps[2]),
            )

        half = 0.5 * self.cell_size * self.repeat_cell
        self._finalize_frep(
            _raw_field,
            (-half[0], half[0], -half[1], half[1], -half[2], half[2]),
        )

    @property
    def raw_field(self: Tpms) -> Field:
        """The raw (non-SDF-normalized) TPMS field callable."""
        return self._raw_field_func

    def as_sheet(self: Tpms, thickness: float | None = None) -> Shape:
        """Return an F-rep Shape representing a TPMS sheet of given thickness.

        Uses the SDF-normalized field, so *thickness* is in physical units.
        If *thickness* is ``None``, uses ``self.offset`` (which may be a
        scalar, an array sampled on ``self.grid``, or a callable in which
        case the callable form is used directly).
        """
        from .implicit_ops import shell  # noqa: PLC0415
        from .shape import Shape  # noqa: PLC0415

        if thickness is not None:
            t: float | npt.NDArray | Field = float(thickness)
        elif self._offset_func is not None:
            t = self._offset_func
        else:
            t = self._offset
        return shell(Shape(func=self._func, bounds=self._bounds), t)

    def _half_offset_field(self: Tpms) -> Field | float:
        """Return half the offset as a callable (variable) or scalar (constant).

        Variable offset stored as a callable can be re-evaluated on the
        marching-cubes grid; an array offset (sampled on ``self.grid``)
        cannot be remapped to arbitrary points, so the only safe fallback
        for arrays is to use a scalar 0 (sheet/skeletal degenerate to the
        zero-isosurface).
        """
        if self._offset_func is not None:
            f = self._offset_func

            def _half(x, y, z, _f=f):  # noqa: ANN001
                return 0.5 * _f(x, y, z)

            return _half
        if isinstance(self._offset, (int, float)):
            return 0.5 * float(self._offset)
        # array — no safe re-evaluation; degenerate to zero (skeletal at f=0).
        return 0.0

    def as_upper_skeletal(self: Tpms) -> Shape:
        """F-rep Shape for the *upper* skeletal: ``{p : f(p) > offset/2}``.

        Volume scales with the chosen offset (smaller offset ⇒ larger
        skeletal), matching the historical CadQuery behaviour and the VTK
        grid-clip path.
        """
        from .implicit_ops import from_field  # noqa: PLC0415

        f = self._func
        h = self._half_offset_field()
        if callable(h):

            def _upper(x, y, z, _f=f, _h=h):  # noqa: ANN001
                return -_f(x, y, z) + _h(x, y, z)

        else:

            def _upper(x, y, z, _f=f, _h=h):  # noqa: ANN001
                return -_f(x, y, z) + _h

        return from_field(func=_upper, bounds=self._bounds)

    def as_lower_skeletal(self: Tpms) -> Shape:
        """F-rep Shape for the *lower* skeletal: ``{p : f(p) < -offset/2}``."""
        from .implicit_ops import from_field  # noqa: PLC0415

        f = self._func
        h = self._half_offset_field()
        if callable(h):

            def _lower(x, y, z, _f=f, _h=h):  # noqa: ANN001
                return _f(x, y, z) + _h(x, y, z)

        else:

            def _lower(x, y, z, _f=f, _h=h):  # noqa: ANN001
                return _f(x, y, z) + _h

        return from_field(func=_lower, bounds=self._bounds)

    def as_surface(self: Tpms) -> Shape:
        """Return F-rep Shape for the (open) zero-isosurface, no thickness.

        Same field as :meth:`as_lower_skeletal`; meant for ``type_part="surface"``.
        Marching cubes will produce an open shell — there is no enclosed volume.
        """
        from .implicit_ops import from_field  # noqa: PLC0415

        return from_field(func=self._func, bounds=self._bounds)

    def _cell_box(self: Tpms) -> Shape:
        """SDF Shape of this TPMS' cell (cell_size × repeat_cell, centered origin)."""
        from .implicit_ops import box  # noqa: PLC0415

        dims = tuple(float(d) for d in (self.cell_size * self.repeat_cell))
        return box(dims=dims, center=(0.0, 0.0, 0.0))

    def _clipped_sheet(self: Tpms) -> Shape:
        """Sheet F-rep clipped to the cell box.

        When the offset approaches its upper limit the sheet field is
        uniformly negative inside the cell and marching cubes finds no
        boundary — clipping by the cell box makes the box face the closed
        boundary, yielding a full-cell mesh as expected.
        """
        from .implicit_ops import intersection  # noqa: PLC0415

        return intersection(self.as_sheet(), self._cell_box())

    def _clipped_upper_skeletal(self: Tpms) -> Shape:
        """Upper skeletal F-rep clipped to the cell box (closed under marching cubes)."""
        from .implicit_ops import intersection  # noqa: PLC0415

        return intersection(self.as_upper_skeletal(), self._cell_box())

    def _clipped_lower_skeletal(self: Tpms) -> Shape:
        """Lower skeletal F-rep clipped to the cell box (closed under marching cubes)."""
        from .implicit_ops import intersection  # noqa: PLC0415

        return intersection(self.as_lower_skeletal(), self._cell_box())

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
        # Reset cached callable form on every assignment.
        self._offset_func = None
        if isinstance(offset, (int, float, np.ndarray)):
            self._offset = offset
        elif isinstance(offset, OffsetGrading):
            self._offset = offset.compute_offset(self.grid)
        elif callable(offset):
            # Keep the callable so the F-rep path can re-evaluate it on the
            # marching-cubes grid (which differs from ``self.grid``).
            self._offset_func = offset
            self._offset = offset(self.grid.x, self.grid.y, self.grid.z).ravel("F")
        else:
            err_msg = "offset must be a float, a numpy array or a callable"
            raise TypeError(err_msg)

        self._update_grid_offset()
        self.offset_updated = True

    def _mesh_to_shell(self: Tpms, mesh: pv.PolyData) -> cq.Shape:
        """Convert a triangulated PyVista mesh to a CadQuery ``Shape``.

        Builds one ``cq.Face`` per triangle, then *sews* them with OCCT's
        ``BRepBuilderAPI_Sewing`` so adjacent triangles share edges — without
        sewing, ``Shell.makeShell`` returns a disjoint compound and the result
        cannot be turned into a Solid.  Returns the sewn shape (a
        ``TopoDS_Shell`` if sewing succeeded into one shell, otherwise the
        compound that sewing produced).
        """
        from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing  # noqa: PLC0415

        if not mesh.is_all_triangles:
            mesh.triangulate(inplace=True)
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

        if not faces:
            err_msg = "Mesh has no triangles to convert into a shell"
            raise ShellCreationError(err_msg)

        sewing = BRepBuilderAPI_Sewing()
        for f in faces:
            sewing.Add(f.wrapped)
        sewing.Perform()
        return cq.Shape(sewing.SewedShape())

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

    _VALID_PARTS = ("sheet", "lower skeletal", "upper skeletal", "surface")

    def _frep_part(self: Tpms, type_part: TpmsPartType) -> Shape:
        """Pick the F-rep :class:`Shape` for *type_part*.

        Skeletals are intersected with the cell box so marching cubes
        produces a *closed* shell (the unclipped skeletal field is
        unbounded).  ``"surface"`` and ``"sheet"`` are already bounded by
        construction (zero-isosurface and `shell()`-clipped, respectively).
        """
        if type_part == "sheet":
            return self._clipped_sheet()
        if type_part == "upper skeletal":
            return self._clipped_upper_skeletal()
        if type_part == "lower skeletal":
            return self._clipped_lower_skeletal()
        if type_part == "surface":
            return self.as_surface()
        err_msg = f"type_part {type_part!r} must be one of {self._VALID_PARTS}"
        raise ValueError(err_msg)

    def _isotropic_resolution(self: Tpms) -> int:
        """Map ``self.resolution`` (per-axis) to an isotropic Shape resolution.

        ``Shape.generate_vtk`` takes a single resolution; we use the geometric
        mean of the per-axis cell counts so total grid points stay proportional.
        """
        return max(int(self.resolution * np.cbrt(np.prod(self.repeat_cell))), 10)

    def generate(
        self: Tpms,
        type_part: TpmsPartType = "sheet",
        smoothing: int = 0,
        algo_resolution: int | None = None,
        **_: KwargsGenerateType,
    ) -> cq.Shape:
        """Generate the OCCT/CadQuery shape of the requested TPMS part.

        Pure F-rep pipeline: pick the SDF Shape via :meth:`_frep_part`, run
        marching cubes through :meth:`Shape.generate_vtk`, optionally smooth,
        then build a CadQuery ``Shell``.  The same SDF + same marching-cubes
        grid is used by :meth:`generate_vtk`, so volumes converge to identical
        values up to discretization.

        :param type_part: ``"sheet"``, ``"lower skeletal"``, ``"upper skeletal"``
            or ``"surface"`` (open zero-isosurface, no thickness)
        :param smoothing: number of Laplacian smoothing iterations on the mesh
        :param algo_resolution: temporary-TPMS resolution for density→offset
            search (only used when ``self.density`` is set)
        """
        if type_part not in self._VALID_PARTS:
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
        else:
            if self.density is not None:
                self._compute_offset_to_fit_density(
                    part_type=type_part,
                    resolution=algo_resolution,
                )
            self._check_offset(type_part)

        frep = self._frep_part(type_part)
        mesh = frep.generate_vtk(
            bounds=self._bounds,
            resolution=self._isotropic_resolution(),
        )
        if smoothing > 0:
            mesh.smooth(n_iter=smoothing, feature_smoothing=True, inplace=True)
            mesh.clean(inplace=True)

        if mesh.n_cells == 0:
            err_msg = (
                f"Marching cubes produced an empty mesh for '{type_part}'; "
                "check offset / density / resolution."
            )
            raise ShellCreationError(err_msg)

        shape = self._mesh_to_shell(mesh)

        if type_part != "surface":
            # Closed shell ⇒ try to upgrade into a Solid so volume queries and
            # downstream booleans behave correctly.  If sewing produced a
            # compound (non-manifold patches) or OCCT refuses, keep the shell.
            shape = self._try_make_solid(shape)

        shape = rotate(obj=shape, center=(0, 0, 0), rotation=self.orientation)
        return shape.translate(self.center)

    @staticmethod
    def _try_make_solid(shape: cq.Shape) -> cq.Shape:
        """Best-effort upgrade of a sewn shell into a closed Solid.

        Returns the original shape unchanged if the sewn result is a Compound
        (multiple disjoint shells, can't be a single solid) or if OCCT refuses
        the conversion.
        """
        from OCP.TopAbs import TopAbs_SHELL  # noqa: PLC0415
        from OCP.TopExp import TopExp_Explorer  # noqa: PLC0415
        from OCP.TopoDS import TopoDS  # noqa: PLC0415

        wrapped = shape.wrapped
        # Already a Shell? Try to make a Solid directly.
        if wrapped.ShapeType() == TopAbs_SHELL:
            try:
                shell = TopoDS.Shell_s(wrapped)
                return cq.Shape(cq.Solid.makeSolid(cq.Shell(shell)).wrapped)
            except (ValueError, RuntimeError):
                return shape

        # Compound: extract Shells, build a Solid per closed shell, fuse.
        exp = TopExp_Explorer(wrapped, TopAbs_SHELL)
        solids: list[cq.Shape] = []
        while exp.More():
            try:
                shell = TopoDS.Shell_s(exp.Current())
                solid = cq.Solid.makeSolid(cq.Shell(shell))
                solids.append(cq.Shape(solid.wrapped))
            except (ValueError, RuntimeError):
                pass
            exp.Next()

        if not solids:
            return shape
        if len(solids) == 1:
            return solids[0]
        return fuseShapes(cqShapeList=solids, retain_edges=False)

    def generate_vtk(
        self: Tpms,
        type_part: TpmsPartType = "sheet",
        algo_resolution: int | None = None,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate the PyVista mesh of the requested TPMS part.

        Same F-rep pipeline as :meth:`generate` (skeletals are clipped to the
        cell box), so the two outputs share the exact same triangulation and
        therefore the same volume.

        :param type_part: ``"sheet"``, ``"lower skeletal"``, ``"upper skeletal"``
            or ``"surface"``
        :param algo_resolution: temporary-TPMS resolution for density→offset
            search (only used when ``self.density`` is set)
        """
        if type_part == "surface":
            return self.surface
        if type_part not in ["sheet", "lower skeletal", "upper skeletal"]:
            err_msg = (
                f"type_part ({type_part}) must be 'sheet', 'lower skeletal', "
                "'upper skeletal' or 'surface'"
            )
            raise ValueError(err_msg)

        if self.density == 1.0:
            # Short-circuit: at full density the part fills the whole cell,
            # return the cell-box mesh directly (exact volume — the
            # marching-cubes path would give a slight discretization gap).
            box = pv.Box(bounds=self._bounds, level=0, quads=False)
            box = box.extract_surface().clean().triangulate()
            box = rotate(box, center=(0, 0, 0), rotation=self.orientation)
            return box.translate(xyz=self.center)

        if self.density is not None:
            self._compute_offset_to_fit_density(
                part_type=type_part,
                resolution=algo_resolution,
            )
        # Note: no `_check_offset` here — the F-rep VTK path handles negative
        # / zero / variable offsets gracefully.  ``generate()`` still applies
        # the historical CAD-side restriction.

        frep = self._frep_part(type_part)
        polydata = frep.generate_vtk(
            bounds=self._bounds,
            resolution=self._isotropic_resolution(),
        )

        polydata = rotate(polydata, center=(0, 0, 0), rotation=self.orientation)
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

    def _setup_frep_field(self: CylindricalTpms) -> None:
        """Build F-rep field in Cartesian coordinates (inverse cylindrical mapping)."""
        k_x, k_y, k_z = 2.0 * np.pi / self.cell_size
        ps = self.phase_shift
        cyl_r = self.cylinder_radius
        unit_theta = self.unit_theta

        def _raw_field(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            rho = np.sqrt(x**2 + y**2) - cyl_r
            theta = np.arctan2(y, x) / unit_theta
            return self.surface_function(
                k_x * (rho + ps[0]),
                k_y * (theta + ps[1]),
                k_z * (z + ps[2]),
            )

        r_max = cyl_r + 0.5 * self.cell_size[0] * self.repeat_cell[0]
        half_z = 0.5 * self.cell_size[2] * self.repeat_cell[2]
        self._finalize_frep(
            _raw_field,
            (-r_max, r_max, -r_max, r_max, -half_z, half_z),
        )


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

    def _setup_frep_field(self: SphericalTpms) -> None:
        """Build F-rep field in Cartesian coordinates (inverse spherical mapping)."""
        k_x, k_y, k_z = 2.0 * np.pi / self.cell_size
        ps = self.phase_shift
        sph_r = self.sphere_radius
        unit_theta = self.unit_theta
        unit_phi = self.unit_phi

        def _raw_field(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            rho_cart = np.sqrt(x**2 + y**2 + z**2)
            rho = rho_cart - sph_r
            theta = (
                np.arccos(np.clip(z / np.maximum(rho_cart, 1e-30), -1, 1))
                - np.pi / 2.0
            ) / unit_theta
            phi = np.arctan2(y, x) / unit_phi
            return self.surface_function(
                k_x * (rho + ps[0]),
                k_y * (theta + ps[1]),
                k_z * (phi + ps[2]),
            )

        r_max = sph_r + 0.5 * self.cell_size[0] * self.repeat_cell[0]
        self._finalize_frep(
            _raw_field,
            (-r_max, r_max, -r_max, r_max, -r_max, r_max),
        )


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

    def _density_envelope_volume(self: Infill) -> float:
        """Density is measured against the input object's volume."""
        return abs(self.obj.volume)

    def _cell_box(self: Infill) -> Shape:
        """Override the cell-box clip with the *object envelope* SDF.

        Without this override the F-rep marching cubes path (used by
        :meth:`Tpms.generate_vtk` and :meth:`Tpms.generate`) would clip the
        TPMS to the cartesian bounding box rather than to the input object —
        producing volumes well above ``obj.volume`` and corrupting density.
        """
        from .implicit_ops import from_field  # noqa: PLC0415

        bounds_arr = np.array(self.obj.bounds)
        bounds: BoundsType = tuple(bounds_arr.tolist())

        envelope = self.obj

        def _envelope_sdf(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            pts = pv.PolyData(np.column_stack([x.ravel(), y.ravel(), z.ravel()]))
            pts.compute_implicit_distance(envelope, inplace=True)
            return np.asarray(pts["implicit_distance"]).reshape(x.shape)

        return from_field(_envelope_sdf, bounds=bounds)


class Conformal(Tpms):
    """TPMS that conforms to the local frame of an envelope surface.

    Where :class:`CylindricalTpms` and :class:`SphericalTpms` use a closed-form
    inverse coordinate map to wrap the TPMS around a known surface (cylinder /
    sphere), :class:`Conformal` does the same for *any* envelope surface
    (PolyData) by using the signed distance to the envelope as the radial
    coordinate.

    Local frame at every Cartesian point ``p``:

    - ``w(p) = signed_distance(p, envelope)`` — *radial* coordinate; ``w=0`` on
      the envelope surface, ``w<0`` inside, ``w>0`` outside.  Computed by
      :meth:`pyvista.PolyData.compute_implicit_distance`.
    - ``u(p) = (p - center) · t``, ``v(p) = (p - center) · b`` — *tangential*
      coordinates expressed in the global frame ``(t, b, n=default_axis)``,
      where ``t`` and ``b`` are the two axes of the plane perpendicular to
      ``default_tangent_axis``.  Default tangent axis is ``(1, 0, 0)``.

    The TPMS field is then ``surface_function(k_w · w, k_u · u, k_v · v)``,
    so unit cells stack along the envelope normal direction (concentric
    "shells" at distance ``±cell_size``, ``±2·cell_size``, …) and tile the
    perpendicular plane in the other two axes.  Compared to a Cartesian
    :class:`Infill`, this produces shells that *follow* the envelope shape
    rather than being sliced by it.

    Tangentially-intrinsic (true conformal) parametrizations require
    per-point local-tangent rotation; that's a follow-up — this class is the
    minimum useful step.
    """

    def __init__(  # noqa: PLR0913
        self: Conformal,
        envelope: pv.PolyData,
        surface_function: Field,
        offset: float | OffsetGrading | Field | None = None,
        cell_size: float | Sequence[float] | npt.NDArray[np.float64] | None = None,
        repeat_cell: int | Sequence[int] | npt.NDArray[np.int8] | None = None,
        phase_shift: Sequence[float] = (0.0, 0.0, 0.0),
        resolution: int = 20,
        density: float | None = None,
        default_tangent_axis: Sequence[float] = (1.0, 0.0, 0.0),
        clip_to_envelope: bool = True,
    ) -> None:
        """Build a TPMS that wraps an envelope surface.

        :param envelope: the surface (``pv.PolyData``) the TPMS will follow.
            Volumes are computed relative to the envelope volume when
            ``clip_to_envelope=True``.
        :param surface_function: tpms function or custom function ``f(x,y,z)=0``
        :param offset: offset / thickness of the TPMS sheet
        :param cell_size: unit cell size; auto-derived from envelope bounds if None
        :param repeat_cell: number of cells per axis; auto-derived if None
        :param phase_shift: phase shift in the (w, u, v) local frame
        :param resolution: per-axis grid resolution
        :param density: target density relative to envelope volume (mutex with
            ``offset``)
        :param default_tangent_axis: world axis used to define the tangential
            plane for ``(u, v)``.  Default ``(1, 0, 0)`` ⇒ ``u = y - cy``,
            ``v = z - cz``.
        :param clip_to_envelope: if ``True``, the resulting grid / mesh is
            clipped by the envelope surface (TPMS only fills the interior).
        """
        self.envelope = envelope
        self.clip_to_envelope = clip_to_envelope

        # Frame: n = default_tangent_axis (radial-projection axis), then build
        # an orthonormal (t, b) basis for the perpendicular plane.
        n = np.asarray(default_tangent_axis, dtype=np.float64)
        n_norm = np.linalg.norm(n)
        if n_norm < 1e-12:
            err_msg = "default_tangent_axis must be a non-zero vector"
            raise ValueError(err_msg)
        n /= n_norm
        # Pick a stable secondary direction not parallel to n.
        helper = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        t = helper - np.dot(helper, n) * n
        t /= np.linalg.norm(t)
        b = np.cross(n, t)
        self._frame_n = n
        self._frame_t = t
        self._frame_b = b
        self.default_tangent_axis = tuple(n.tolist())

        # Auto-derive cell_size / repeat_cell from envelope bbox if missing.
        bounds = np.array(envelope.bounds)
        margin = 1.001
        envelope_dim = margin * (bounds[1::2] - bounds[::2])

        if cell_size is not None and repeat_cell is not None:
            err_msg = (
                "cell_size and repeat_cell cannot be given at the same time, "
                "one is computed from the other."
            )
            raise ValueError(err_msg)
        if cell_size is None and repeat_cell is None:
            repeat_cell = (4, 4, 4)
        if cell_size is not None:
            repeat_cell = np.maximum(
                np.round(envelope_dim / np.asarray(cell_size, dtype=float)).astype(int),
                1,
            )
        elif repeat_cell is not None:
            cell_size = envelope_dim / np.asarray(repeat_cell, dtype=float)

        self._envelope_center = np.asarray(envelope.center, dtype=np.float64)

        super().__init__(
            surface_function=surface_function,
            offset=offset,
            phase_shift=phase_shift,
            cell_size=cell_size,
            repeat_cell=repeat_cell,
            resolution=resolution,
            density=density,
        )

    # -- Local-frame field -------------------------------------------------

    def _envelope_distance(
        self: Conformal,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Vectorized signed distance to the envelope (negative inside)."""
        pts = pv.PolyData(np.column_stack([x.ravel(), y.ravel(), z.ravel()]))
        pts.compute_implicit_distance(self.envelope, inplace=True)
        return np.asarray(pts["implicit_distance"]).reshape(x.shape)

    def _local_coords(
        self: Conformal,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> tuple[
        npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]
    ]:
        """Project Cartesian ``(x, y, z)`` onto the conformal local frame.

        - ``w`` = signed distance to the envelope (radial coord)
        - ``u`` = ``(p - envelope_center) · t``
        - ``v`` = ``(p - envelope_center) · b``
        """
        cx, cy, cz = self._envelope_center
        dx, dy, dz = x - cx, y - cy, z - cz
        u = dx * self._frame_t[0] + dy * self._frame_t[1] + dz * self._frame_t[2]
        v = dx * self._frame_b[0] + dy * self._frame_b[1] + dz * self._frame_b[2]
        w = self._envelope_distance(x, y, z)
        return u, v, w

    def _setup_frep_field(self: Conformal) -> None:
        """Build F-rep field that evaluates the surface function in (w, u, v)."""
        k_x, k_y, k_z = 2.0 * np.pi / self.cell_size
        ps = self.phase_shift

        def _raw_field(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            u, v, w = self._local_coords(x, y, z)
            return self.surface_function(
                k_x * (w + ps[0]),
                k_y * (u + ps[1]),
                k_z * (v + ps[2]),
            )

        bounds_arr = np.array(self.envelope.bounds)
        bounds: BoundsType = tuple(bounds_arr.tolist())
        self._finalize_frep(_raw_field, bounds)

    # -- Override grid creation to clip to the envelope --------------------

    def _create_grid(
        self: Conformal,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> pv.StructuredGrid:
        """Cartesian grid covering the envelope bbox, optionally clipped."""
        cx, cy, cz = self._envelope_center
        grid = super()._create_grid(x + cx, y + cy, z + cz)
        if self.clip_to_envelope:
            grid = grid.clip_surface(self.envelope)
            logging.info(
                "Conformal grid clipped to envelope: %d points", grid.n_points,
            )
        return grid

    def _density_envelope_volume(self: Conformal) -> float:
        """Density is measured against the envelope volume."""
        return abs(self.envelope.volume)

    # -- Override the F-rep "cell box" with the envelope itself ------------

    def _cell_box(self: Conformal) -> Shape:
        """SDF Shape of the envelope (used to clip TPMS parts to the interior).

        Replaces the parent's axis-aligned cell-box clip — the natural cell of
        a Conformal TPMS *is* the envelope.  Returns a Shape whose field is
        ``signed_distance(p, envelope)`` (negative inside).
        """
        from .implicit_ops import from_field  # noqa: PLC0415

        bounds_arr = np.array(self.envelope.bounds)
        bounds: BoundsType = tuple(bounds_arr.tolist())

        def _envelope_sdf(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            return self._envelope_distance(x, y, z)

        return from_field(_envelope_sdf, bounds=bounds)


# Re-export for backward compatibility
from .shape import BoundsType, ShellCreationError  # noqa: E402

__all__ = [
    "Conformal",
    "CylindricalTpms",
    "Infill",
    "ShellCreationError",
    "SphericalTpms",
    "Tpms",
]
