"""
TPMS.

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
from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING, Literal

import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.optimize import root_scalar

from microgen.operations import fuse_shapes, rotate

from .shape import Shape

if TYPE_CHECKING:
    from microgen.cad import CadShape
    from microgen.shape import KwargsGenerateType, TpmsPartType, Vector3DType
from .tpms_grading import OffsetGrading

logging.basicConfig(level=logging.INFO)
Field = Callable[
    [npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]],
    npt.NDArray[np.float64],
]

_DIM = 3


def _wedge_sdf_2d(
    x: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
    half_angle: float,
) -> npt.NDArray[np.float64]:
    """
    SDF of a 2D wedge centred on +X with half-aperture ``half_angle``.

    Negative inside, positive outside, zero on the bounding rays.
    Z-independent — used as a clipping primitive for cylindrical / spherical
    partial wraps.
    """
    return np.cos(half_angle) * np.abs(y) - np.sin(half_angle) * x


def _double_cone_sdf(
    x: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
    z: npt.NDArray[np.float64],
    half_polar: float,
) -> npt.NDArray[np.float64]:
    """
    SDF of a double cone around the Z axis with half-aperture ``half_polar``.

    Symmetric about the XY plane; points inside iff their angle from ±Z is
    less than ``half_polar``.  Used to clip spherical partial-θ coverage.
    """
    rho_xy = np.sqrt(x * x + y * y)
    return np.cos(half_polar) * rho_xy - np.sin(half_polar) * np.abs(z)


def _interp_along_curve(
    query_s: npt.NDArray[np.float64],
    curve_s: npt.NDArray[np.float64],
    values: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """
    Per-axis ``np.interp`` of an ``(M, 3)`` curve attribute at ``query_s``.

    Returns an ``(N, 3)`` array.  Wraps the three-axis loop so callers don't
    repeat the pattern.
    """
    return np.column_stack(
        [np.interp(query_s, curve_s, values[:, k]) for k in range(values.shape[1])],
    )


class Tpms(Shape):
    """
    Triply Periodical Minimal Surfaces.

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
        r"""
        Class used to generate TPMS geometries (sheet or skeletals parts).

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
            offset
            if (
                offset is not None
                and callable(offset)
                and not isinstance(offset, OffsetGrading)
            )
            else None
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
        """
        Return the offset corresponding to the required density.

        :param surface_function: tpms function
        :param part_type: type of the part (sheet, lower skeletal or upper skeletal)
        :param density: Required density, 0.5 for 50%
        :param resolution: resolution of the tpms used to compute the offset

        :return: corresponding offset value
        """
        return Tpms(
            surface_function=surface_function,
            density=density,
        )._compute_offset_to_fit_density(part_type=part_type, resolution=resolution)

    def _density_envelope_volume(self: Tpms) -> float:
        """
        Return the *envelope* volume that ``density`` is measured against.

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
        """
        Compute the offset that yields the requested density.

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
                strict=False,
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
        from .implicit_ops import from_field, normalize_to_sdf

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
        """
        Return an F-rep Shape representing a TPMS sheet of given thickness.

        Uses the SDF-normalized field, so *thickness* is in physical units.
        If *thickness* is ``None``, uses ``self.offset`` (which may be a
        scalar, an array sampled on ``self.grid``, or a callable in which
        case the callable form is used directly).
        """
        from .implicit_ops import shell
        from .shape import Shape

        if thickness is not None:
            t: float | npt.NDArray | Field = float(thickness)
        elif self._offset_func is not None:
            t = self._offset_func
        else:
            t = self._offset
        return shell(Shape(func=self._func, bounds=self._bounds), t)

    def _half_offset_field(self: Tpms) -> Field | float:
        """
        Return half the offset as a callable (variable) or scalar (constant).

        Variable offset stored as a callable can be re-evaluated on the
        marching-cubes grid; an array offset (sampled on ``self.grid``)
        cannot be remapped to arbitrary points, so the only safe fallback
        for arrays is to use a scalar 0 (sheet/skeletal degenerate to the
        zero-isosurface).
        """
        if self._offset_func is not None:
            f = self._offset_func

            def _half(x, y, z, _f=f):
                return 0.5 * _f(x, y, z)

            return _half
        if isinstance(self._offset, (int, float)):
            return 0.5 * float(self._offset)
        # array — no safe re-evaluation; degenerate to zero (skeletal at f=0).
        return 0.0

    def as_upper_skeletal(self: Tpms) -> Shape:
        """
        F-rep Shape for the *upper* skeletal: ``{p : f(p) > offset/2}``.

        Volume scales with the chosen offset (smaller offset ⇒ larger
        skeletal), matching the historical CadQuery behaviour and the VTK
        grid-clip path.
        """
        from .implicit_ops import from_field

        f = self._func
        h = self._half_offset_field()
        if callable(h):

            def _upper(x, y, z, _f=f, _h=h):
                return -_f(x, y, z) + _h(x, y, z)

        else:

            def _upper(x, y, z, _f=f, _h=h):
                return -_f(x, y, z) + _h

        return from_field(func=_upper, bounds=self._bounds)

    def as_lower_skeletal(self: Tpms) -> Shape:
        """F-rep Shape for the *lower* skeletal: ``{p : f(p) < -offset/2}``."""
        from .implicit_ops import from_field

        f = self._func
        h = self._half_offset_field()
        if callable(h):

            def _lower(x, y, z, _f=f, _h=h):
                return _f(x, y, z) + _h(x, y, z)

        else:

            def _lower(x, y, z, _f=f, _h=h):
                return _f(x, y, z) + _h

        return from_field(func=_lower, bounds=self._bounds)

    def as_surface(self: Tpms) -> Shape:
        """
        Return F-rep Shape for the (open) zero-isosurface, no thickness.

        Same field as :meth:`as_lower_skeletal`; meant for ``type_part="surface"``.
        Marching cubes will produce an open shell — there is no enclosed volume.
        """
        from .implicit_ops import from_field

        return from_field(func=self._func, bounds=self._bounds)

    def _cell_box(self: Tpms) -> Shape:
        """SDF Shape of this TPMS' cell (cell_size × repeat_cell, centered origin)."""
        from .implicit_ops import box

        dims = tuple(float(d) for d in (self.cell_size * self.repeat_cell))
        return box(dims=dims, center=(0.0, 0.0, 0.0))

    def _clipped_sheet(self: Tpms) -> Shape:
        """
        Sheet F-rep clipped to the cell box.

        When the offset approaches its upper limit the sheet field is
        uniformly negative inside the cell and marching cubes finds no
        boundary — clipping by the cell box makes the box face the closed
        boundary, yielding a full-cell mesh as expected.
        """
        from .implicit_ops import intersection

        return intersection(self.as_sheet(), self._cell_box())

    def _clipped_upper_skeletal(self: Tpms) -> Shape:
        """Upper skeletal F-rep clipped to the cell box (closed under marching cubes)."""
        from .implicit_ops import intersection

        return intersection(self.as_upper_skeletal(), self._cell_box())

    def _clipped_lower_skeletal(self: Tpms) -> Shape:
        """Lower skeletal F-rep clipped to the cell box (closed under marching cubes)."""
        from .implicit_ops import intersection

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
            # Sample the callable on ``self.grid`` for legacy / property-path
            # consumers (``grid_sheet`` etc.).  StructuredGrid exposes
            # ``.x/.y/.z`` meshgrids; UnstructuredGrid (e.g. after
            # ``clip_surface`` in :class:`Infill`) does not — use raw points.
            if hasattr(self.grid, "x"):
                vals = offset(self.grid.x, self.grid.y, self.grid.z).ravel("F")
            else:
                pts = np.asarray(self.grid.points)
                vals = offset(pts[:, 0], pts[:, 1], pts[:, 2])
            self._offset = vals
        else:
            err_msg = "offset must be a float, a numpy array or a callable"
            raise TypeError(err_msg)

        self._update_grid_offset()
        self.offset_updated = True

    def _mesh_to_shell(self: Tpms, mesh: pv.PolyData) -> CadShape:
        """Convert a triangulated PyVista mesh to an OCCT ``CadShape``.

        Delegates to :func:`microgen.cad.mesh_to_shell_brep` (one planar BREP
        face per triangle), the slower-but-correct path that produces a shell
        with real surface geometry — required because downstream callers may
        use the result as a cutting tool in boolean ops.
        """
        from microgen.cad import ShellCreationError as _CadShellError  # noqa: PLC0415
        from microgen.cad import mesh_to_shell_brep  # noqa: PLC0415

        if not mesh.is_all_triangles:
            mesh.triangulate(inplace=True)
        triangles = mesh.faces.reshape(-1, 4)[:, 1:]
        points = np.asarray(mesh.points, dtype=np.float64)

        try:
            return mesh_to_shell_brep(points, triangles)
        except _CadShellError as err:
            err_msg = (
                "Failed to create the shell; "
                "try to increase the resolution or the smoothing."
            )
            raise ShellCreationError(err_msg) from err

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
        """
        Pick the F-rep :class:`Shape` for *type_part*.

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
        """
        Map ``self.resolution`` (per-axis) to an isotropic Shape resolution.

        ``Shape.generate_vtk`` takes a single resolution; we use the geometric
        mean of the per-axis cell counts so total grid points stay proportional.
        """
        return max(int(self.resolution * np.cbrt(np.prod(self.repeat_cell))), 10)

    def _envelope_mesh_at_full_density(self: Tpms) -> pv.PolyData:
        """
        Return the cell-envelope mesh used by the ``density=1.0`` shortcut.

        Plain :class:`Tpms` has an axis-aligned box envelope, returned exactly
        via :class:`pyvista.Box`.  Subclasses with non-cubic envelopes
        (cylindrical / spherical shells, ``Infill`` object, etc.) override
        this to fall back to marching cubes on their ``_cell_box`` SDF — see
        :meth:`_envelope_mesh_via_cell_box`.
        """
        return (
            pv.Box(bounds=self._bounds, level=0, quads=False)
            .extract_surface()
            .clean()
            .triangulate()
        )

    def _envelope_mesh_via_cell_box(self: Tpms) -> pv.PolyData:
        """
        Marching-cubes mesh of the (potentially non-cubic) ``_cell_box``.

        Used by subclasses that override :meth:`_envelope_mesh_at_full_density`
        to delegate to the F-rep envelope SDF.
        """
        envelope_shape = self._cell_box()
        return envelope_shape.generate_vtk(
            bounds=envelope_shape.bounds or self._bounds,
            resolution=self._isotropic_resolution(),
        )

    def generate(
        self: Tpms,
        type_part: TpmsPartType = "sheet",
        smoothing: int = 0,
        algo_resolution: int | None = None,
        **_: KwargsGenerateType,
    ) -> CadShape:
        """Generate an OCCT CAD shape of the requested TPMS part.

        Pure F-rep pipeline: pick the SDF Shape via :meth:`_frep_part`, run
        marching cubes through :meth:`Shape.generate_vtk`, optionally smooth,
        then build an OCCT ``Shell`` via :func:`microgen.cad.mesh_to_shell_brep`.
        The same SDF + same marching-cubes grid is used by :meth:`generate_vtk`,
        so volumes converge to identical values up to discretisation.

        Requires the optional ``[cad]`` install extra (``cadquery-ocp``).

        :param type_part: ``"sheet"``, ``"lower skeletal"``, ``"upper skeletal"``
            or ``"surface"`` (open zero-isosurface, no thickness)
        :param smoothing: number of Laplacian smoothing iterations on the mesh
        :param algo_resolution: temporary-TPMS resolution for density→offset
            search (only used when ``self.density`` is set)
        :return: :class:`microgen.cad.CadShape` wrapping an OCCT ``TopoDS_Shell``
            (or a ``TopoDS_Solid`` when sewing succeeded into a closed shell).
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
    def _try_make_solid(shape: CadShape) -> CadShape:
        """Best-effort upgrade of a sewn shell into a closed Solid.

        Returns the original shape unchanged if the sewn result is a Compound
        (multiple disjoint shells, can't be a single solid) or if OCCT refuses
        the conversion.
        """
        from microgen.cad import CadShape as _CadShape  # noqa: PLC0415
        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeSolid  # noqa: PLC0415
        from OCP.TopAbs import TopAbs_SHELL  # noqa: PLC0415
        from OCP.TopExp import TopExp_Explorer  # noqa: PLC0415
        from OCP.TopoDS import TopoDS  # noqa: PLC0415

        wrapped = shape.wrapped
        # Already a Shell? Try to make a Solid directly.
        if wrapped.ShapeType() == TopAbs_SHELL:
            try:
                shell = TopoDS.Shell_s(wrapped)
                return _CadShape(BRepBuilderAPI_MakeSolid(shell).Solid())
            except (ValueError, RuntimeError):
                return shape

        # Compound: extract Shells, build a Solid per closed shell, fuse.
        exp = TopExp_Explorer(wrapped, TopAbs_SHELL)
        solids: list[CadShape] = []
        while exp.More():
            try:
                shell = TopoDS.Shell_s(exp.Current())
                solids.append(_CadShape(BRepBuilderAPI_MakeSolid(shell).Solid()))
            except (ValueError, RuntimeError):
                pass
            exp.Next()

        if not solids:
            return shape
        if len(solids) == 1:
            return solids[0]
        return fuse_shapes(solids, retain_edges=False)

    def generate_vtk(
        self: Tpms,
        type_part: TpmsPartType = "sheet",
        algo_resolution: int | None = None,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """
        Generate the PyVista mesh of the requested TPMS part.

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
            envelope_mesh = self._envelope_mesh_at_full_density()
            envelope_mesh = rotate(
                envelope_mesh,
                center=(0, 0, 0),
                rotation=self.orientation,
            )
            return envelope_mesh.translate(xyz=self.center)

        if self.density is not None:
            self._compute_offset_to_fit_density(
                part_type=type_part,
                resolution=algo_resolution,
            )
        # Note: no `_check_offset` here — the F-rep VTK path handles negative
        # / zero / variable offsets gracefully.  ``generate()`` still applies
        # the historical CAD-side restriction.

        # Subclasses with a *parametric* coordinate frame (CylindricalTpms,
        # SphericalTpms, Sweep) cannot use the F-rep marching-cubes path —
        # the TPMS field has rapid angular oscillations that an isotropic
        # Cartesian grid samples poorly (visible as dotted holes near the
        # pole / axis).  They opt into the legacy grid-clip path by setting
        # ``_uses_parametric_grid = True``: clip the parametric structured
        # grid by the relevant scalar threshold, then extract the surface.
        if getattr(self, "_uses_parametric_grid", False):
            grid_attr = f"grid_{type_part.replace(' ', '_')}"
            # Merge points within an absolute tolerance large enough to close
            # the angular seam (φ=±π for sphere, θ=±π for cylinder) where
            # the structured grid's two boundary faces coincide in Cartesian
            # space up to ~1e-3 floating-point drift, but small enough not
            # to collapse legitimate cell edges (smallest grid spacing is
            # ~ ``cell_size / resolution``).
            seam_tol = min(
                5e-3,
                0.1 * float(np.min(self.cell_size)) / float(self.resolution),
            )
            polydata = (
                getattr(self, grid_attr)
                .extract_surface()
                .clean(tolerance=seam_tol, absolute=True)
                .triangulate()
            )
        else:
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
        """Deprecated. Use :meth:`generate_vtk` instead."""
        return self.generate_vtk(
            type_part=type_part,
            algo_resolution=algo_resolution,
        )


class CylindricalTpms(Tpms):
    """Class used to generate cylindrical TPMS geometries (sheet or skeletals parts)."""

    # Use the parametric structured-grid clip for ``generate_vtk`` instead of
    # F-rep marching cubes.  An isotropic Cartesian MC grid samples the
    # angular axis poorly near the cylinder axis (rapid θ-derivative ⇒
    # under-sampled holes).  The structured grid in (ρ, θ, z) automatically
    # maps to a Cartesian grid that is dense near the axis.
    _uses_parametric_grid = True

    _envelope_mesh_at_full_density = Tpms._envelope_mesh_via_cell_box

    def __init__(
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
        r"""
        Cylindrical TPMS geometry.

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
        """
        Return the structured cylindrical grid of the TPMS.

        ``y`` carries arc-length along the equator, so the conversion to an
        angle is ``theta = y / radius``.  The earlier ``y * unit_theta`` form
        over-wrapped past ``2π`` by a factor of ``unit_theta * radius``
        (which is only 1 when the user-supplied ``cell_size[1] == 1``),
        leaving a visible wedge of unrendered geometry on one meridian.
        """
        rho = x + self.cylinder_radius
        theta = y / self.cylinder_radius

        grid = pv.StructuredGrid(rho * np.cos(theta), rho * np.sin(theta), z)
        grid["coords"] = np.c_[
            x.ravel(order="F"),
            y.ravel(order="F"),
            z.ravel(order="F"),
        ]
        return grid

    def _shell_extents(self: CylindricalTpms) -> tuple[float, float, float]:
        """Return ``(r_inner, r_outer, half_z)`` of the TPMS-bearing shell."""
        cyl_r = float(self.cylinder_radius)
        delta_r = 0.5 * float(self.cell_size[0]) * float(self.repeat_cell[0])
        r_outer = cyl_r + delta_r
        r_inner = max(cyl_r - delta_r, 0.0)
        half_z = 0.5 * float(self.cell_size[2]) * float(self.repeat_cell[2])
        return r_inner, r_outer, half_z

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

        # Tightened bounds: only ``r_outer`` reach in xy, only ``half_z`` in z.
        # The previous wider Cartesian box wasted ~50 % of MC samples on the
        # interior of the cylinder where the TPMS doesn't live.
        _, r_outer, half_z = self._shell_extents()
        self._finalize_frep(
            _raw_field,
            (-r_outer, r_outer, -r_outer, r_outer, -half_z, half_z),
        )

    def _cell_box(self: CylindricalTpms) -> Shape:
        """
        Cylindrical-shell SDF in Cartesian space (replaces the parent's
        axis-aligned-box SDF, which was geometrically wrong for this class).

        Without this override the F-rep marching cubes path
        (:meth:`Tpms._clipped_sheet`) would clip the TPMS by an axis-aligned
        box of ``cell_size × repeat_cell`` *interpreted as Cartesian* — but
        ``cell_size[1]`` / ``repeat_cell[1]`` are angular (parametric) for
        this class, so the clip discarded ~95 % of the physical shell.

        The shell SDF is::

            radial = sqrt(x² + y²)
            ring   = max(r_inner − radial, radial − r_outer)
            sdf    = max(ring, |z| − half_z)

        Partial wrap (``repeat_cell[1] < n_repeat_to_full_circle``) is
        handled by intersecting the ring with an angular wedge centred on
        the +X axis.
        """
        from .implicit_ops import from_field, intersection

        r_inner, r_outer, half_z = self._shell_extents()
        bounds: BoundsType = (
            -r_outer,
            r_outer,
            -r_outer,
            r_outer,
            -half_z,
            half_z,
        )

        def _shell_sdf(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            radial = np.sqrt(x * x + y * y)
            outer = radial - r_outer
            if r_inner > 0:
                ring = np.maximum(r_inner - radial, outer)
            else:
                ring = outer  # solid disc — no inner wall
            return np.maximum(ring, np.abs(z) - half_z)

        shell = from_field(_shell_sdf, bounds=bounds)

        # Partial angular wrap → intersect with a wedge centred on +X.
        n_full = int(round(2.0 * np.pi / self.unit_theta))
        if int(self.repeat_cell[1]) < n_full:
            half_angle = np.pi * float(self.repeat_cell[1]) / float(n_full)
            wedge = from_field(
                lambda x, y, z, _h=half_angle: _wedge_sdf_2d(x, y, _h),
                bounds=bounds,
            )
            shell = intersection(shell, wedge)

        return shell

    def _isotropic_resolution(self: CylindricalTpms) -> int:
        """
        Per-axis resolution count for the F-rep marching-cubes grid.

        Override the parent's geometric-mean of ``repeat_cell`` (which mixes
        the angular axis with linear ones).  Use only the radial (``[0]``)
        and axial (``[2]``) physical extents, both in world units.
        """
        radial_extent = float(self.cell_size[0]) * float(self.repeat_cell[0])
        axial_extent = float(self.cell_size[2]) * float(self.repeat_cell[2])
        cell_size_min = max(
            min(float(self.cell_size[0]), float(self.cell_size[2])),
            1e-12,
        )
        n = int(self.resolution * max(radial_extent, axial_extent) / cell_size_min)
        return max(n, 10)


class SphericalTpms(Tpms):
    """Class used to generate spherical TPMS geometries (sheet or skeletals parts)."""

    # Same rationale as :class:`CylindricalTpms` — the parametric (ρ, θ, φ)
    # grid stretches naturally in Cartesian space, sampling poles and equator
    # equally well, while uniform Cartesian MC produces pole pathologies.
    _uses_parametric_grid = True

    _envelope_mesh_at_full_density = Tpms._envelope_mesh_via_cell_box

    def __init__(
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
        r"""
        Spherical TPMS geometry.

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
        """
        Return the structured spherical grid of the TPMS.

        ``y`` and ``z`` carry equatorial arc length, so the conversion to
        angles is ``theta = y / radius`` and ``phi = z / radius``.  The
        earlier ``* unit_theta`` / ``* unit_phi`` forms over-wrapped past
        ``π`` / ``2π`` by a factor of ``unit_* * radius`` (which is only
        1 when the user-supplied ``cell_size == 1``), leaving a visible
        wedge of unrendered geometry on one meridian and tiny over-shoots
        at the poles.
        """
        rho = x + self.sphere_radius
        theta = y / self.sphere_radius + np.pi / 2.0
        phi = z / self.sphere_radius

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
                np.arccos(np.clip(z / np.maximum(rho_cart, 1e-30), -1, 1)) - np.pi / 2.0
            ) / unit_theta
            phi = np.arctan2(y, x) / unit_phi
            return self.surface_function(
                k_x * (rho + ps[0]),
                k_y * (theta + ps[1]),
                k_z * (phi + ps[2]),
            )

        # Tightened bounds: only ``r_outer`` cube, not ``r_max + R``.
        _, r_outer = self._shell_radii()
        self._finalize_frep(
            _raw_field,
            (-r_outer, r_outer, -r_outer, r_outer, -r_outer, r_outer),
        )

    def _shell_radii(self: SphericalTpms) -> tuple[float, float]:
        """Return ``(r_inner, r_outer)`` of the spherical TPMS-bearing shell."""
        sph_r = float(self.sphere_radius)
        delta_r = 0.5 * float(self.cell_size[0]) * float(self.repeat_cell[0])
        return max(sph_r - delta_r, 0.0), sph_r + delta_r

    def _cell_box(self: SphericalTpms) -> Shape:
        """
        Spherical-shell SDF in Cartesian space.

        Replaces the parent's axis-aligned-box clip — same fix as
        :meth:`CylindricalTpms._cell_box`, see that docstring for the
        rationale.

        Partial coverage in θ (cone clip on +Z) and φ (wedge clip on +X
        in xy-plane) is supported when ``repeat_cell[1]`` /
        ``repeat_cell[2]`` are smaller than the auto-fill values.
        """
        from .implicit_ops import from_field, intersection

        r_inner, r_outer = self._shell_radii()
        bounds: BoundsType = (
            -r_outer,
            r_outer,
            -r_outer,
            r_outer,
            -r_outer,
            r_outer,
        )

        def _shell_sdf(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            r = np.sqrt(x * x + y * y + z * z)
            outer = r - r_outer
            if r_inner > 0:
                return np.maximum(r_inner - r, outer)
            return outer  # solid ball if inner radius collapsed

        shape = from_field(_shell_sdf, bounds=bounds)

        # Partial θ coverage → double-cone clip around the Z axis.
        n_full_theta = int(round(np.pi / self.unit_theta))
        if int(self.repeat_cell[1]) < n_full_theta:
            half_polar = np.pi * float(self.repeat_cell[1]) / float(n_full_theta)
            shape = intersection(
                shape,
                from_field(
                    lambda x, y, z, _h=half_polar: _double_cone_sdf(x, y, z, _h),
                    bounds=bounds,
                ),
            )

        # Partial φ coverage → wedge clip around the +X axis.
        n_full_phi = int(round(2.0 * np.pi / self.unit_phi))
        if int(self.repeat_cell[2]) < n_full_phi:
            half_azi = np.pi * float(self.repeat_cell[2]) / float(n_full_phi)
            shape = intersection(
                shape,
                from_field(
                    lambda x, y, z, _h=half_azi: _wedge_sdf_2d(x, y, _h),
                    bounds=bounds,
                ),
            )

        return shape

    def _isotropic_resolution(self: SphericalTpms) -> int:
        """
        Use only the radial physical extent (the angular axes are
        parametric and would mislead a geometric mean across them).
        """
        radial_extent = float(self.cell_size[0]) * float(self.repeat_cell[0])
        radius_envelope = self.sphere_radius + radial_extent
        cell_size_radial = max(float(self.cell_size[0]), 1e-12)
        n = int(self.resolution * 2.0 * radius_envelope / cell_size_radial)
        return max(n, 10)


class Sweep(Tpms):
    """
    TPMS along an arbitrary 3D curve — generalisation of CylindricalTpms.

    The TPMS is generated in a tube of radius ``radial_max`` around an
    arbitrary curve, in local coordinates ``(s, r, θ)`` where:

    - ``s`` = arc length along the curve from its start
    - ``r`` = perpendicular distance from the curve at the closest point
    - ``θ`` = angle around the curve in the local tangent-frame plane

    The local frame ``(T, N, B)`` is built once at curve setup time using
    parallel transport along the polyline (tangent from finite differences,
    normal from a seed transported by Rodrigues rotation between adjacent
    samples; closed loops get a holonomy correction distributed linearly).

    The mesh is generated via the **parametric structured-grid path**
    (same as :class:`CylindricalTpms` / :class:`SphericalTpms`): we build a
    structured grid in (s, r, θ) parametric space, evaluate the TPMS field
    on that grid, and map the parametric points to Cartesian using the
    parallel-transport frames.  ``generate_vtk`` then clips the structured
    grid by the relevant scalar threshold — much cleaner than F-rep MC at
    a uniform Cartesian resolution, which would under-sample the angular
    direction and produce dotted artefacts.

    :class:`CylindricalTpms` is a special case of this class (curve = a
    straight line).
    """

    # Use the parametric grid-clip path (same as Cylindrical/Spherical).
    _uses_parametric_grid = True

    _envelope_mesh_at_full_density = Tpms._envelope_mesh_via_cell_box

    def __init__(
        self: Sweep,
        curve_points: npt.NDArray[np.float64]
        | Callable[[float], npt.NDArray[np.float64]],
        surface_function: Field,
        radial_max: float,
        offset: float | OffsetGrading | Field | None = None,
        cell_size: float | Sequence[float] | npt.NDArray[np.float64] = 1.0,
        repeat_cell: int | Sequence[int] | npt.NDArray[np.int8] = 1,
        phase_shift: Sequence[float] = (0.0, 0.0, 0.0),
        resolution: int = 20,
        density: float | None = None,
        seed_normal: Sequence[float] | None = None,
        n_curve_samples: int = 200,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType = (0, 0, 0),
    ) -> None:
        r"""
        Build a TPMS swept along a curve.

        :param curve_points: either an ``(M, 3)`` array of polyline samples
            or a callable ``t \in [0, 1] -> (3,)``.  Callables are sampled
            at ``n_curve_samples`` points before processing.
        :param surface_function: TPMS function ``f(x, y, z)``
        :param radial_max: outer tube radius
        :param offset: TPMS sheet thickness
        :param cell_size: ``(s, r, θ)`` cell size — third axis is angular
            (radians per cell), so a sensible default is to leave it at
            ``1.0`` and tune via ``repeat_cell[2]``.
        :param repeat_cell: ``(n_s, n_r, n_θ)`` — radial cell count is
            usually ``1`` for a thin tube; ``n_θ`` controls how many
            angular cells around the curve.
        :param phase_shift: TPMS phase shift in (s, r, θ)
        :param resolution: per-axis MC grid resolution
        :param density: mutex with ``offset``; density relative to the
            tube volume
        :param seed_normal: initial normal direction for parallel transport;
            defaults to a vector perpendicular to the first tangent
        :param n_curve_samples: number of samples to use when resampling a
            ``Callable`` curve (ignored for polyline input)
        :param center: center of the geometry
        :param orientation: orientation of the geometry
        """
        # Discretise the curve.
        if callable(curve_points):
            ts = np.linspace(0.0, 1.0, int(n_curve_samples))
            curve = np.asarray([curve_points(t) for t in ts], dtype=np.float64)
        else:
            curve = np.asarray(curve_points, dtype=np.float64)
        if curve.ndim != 2 or curve.shape[1] != 3 or curve.shape[0] < 2:
            err_msg = (
                f"curve_points must be an (M, 3) array with M ≥ 2, "
                f"got shape {curve.shape}"
            )
            raise ValueError(err_msg)

        self.curve = curve
        self.radial_max = float(radial_max)

        # Build the parallel-transport frames + arc-length parametrisation.
        self._build_curve_frames(seed_normal=seed_normal)

        # Initialise like a regular TPMS — cell_size now lives in (s, r, θ)
        # parametric space, ``_setup_frep_field`` will use the local frames.
        super().__init__(
            surface_function=surface_function,
            offset=offset,
            phase_shift=phase_shift,
            cell_size=cell_size,
            repeat_cell=repeat_cell,
            resolution=resolution,
            density=density,
            center=center,
            orientation=orientation,
        )

    # -- Curve preprocessing -----------------------------------------------

    def _build_curve_frames(
        self: Sweep,
        seed_normal: Sequence[float] | None,
    ) -> None:
        """
        Pre-compute arc length, tangents, and parallel-transported normals.

        Stores on ``self``:
          ``_curve_s``       arc length per sample, shape (M,)
          ``_curve_T``       unit tangent per sample, shape (M, 3)
          ``_curve_N``       unit normal (parallel-transported), shape (M, 3)
          ``_curve_kdtree``  scipy cKDTree on the curve points
        """
        from scipy.spatial import cKDTree

        pts = self.curve
        m = pts.shape[0]

        # Arc-length parameterisation.
        seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        s_cum = np.concatenate([[0.0], np.cumsum(seg)])
        self._curve_s = s_cum
        self._curve_total_length = float(s_cum[-1])

        # Tangent: forward differences, last point copies its predecessor.
        diffs = np.diff(pts, axis=0)
        diff_norm = np.linalg.norm(diffs, axis=1, keepdims=True)
        diff_norm = np.where(diff_norm > 1e-12, diff_norm, 1.0)
        T = diffs / diff_norm
        T = np.vstack([T, T[-1:]])  # pad final tangent
        self._curve_T = T

        # Seed normal: project the user's seed (or world-up / world-x as
        # fallback) onto the plane perpendicular to T[0].
        if seed_normal is None:
            world_up = np.array([0.0, 0.0, 1.0])
            seed_proj = world_up - np.dot(world_up, T[0]) * T[0]
            if np.linalg.norm(seed_proj) < 1e-6:
                # T[0] is along z — use world-x instead.
                world_x = np.array([1.0, 0.0, 0.0])
                seed_proj = world_x - np.dot(world_x, T[0]) * T[0]
        else:
            seed_proj = np.asarray(seed_normal, dtype=np.float64)
            seed_proj = seed_proj - np.dot(seed_proj, T[0]) * T[0]
            if np.linalg.norm(seed_proj) < 1e-6:
                err_msg = "seed_normal is parallel to the initial tangent"
                raise ValueError(err_msg)
        N = np.empty_like(T)
        N[0] = seed_proj / np.linalg.norm(seed_proj)

        # Parallel transport: rotate previous N around (T_{i-1} × T_i) by
        # the angle between successive tangents (Rodrigues' rotation).
        for i in range(1, m):
            t_prev, t_curr = T[i - 1], T[i]
            cos_angle = float(np.clip(np.dot(t_prev, t_curr), -1.0, 1.0))
            if cos_angle > 1.0 - 1e-9:
                N[i] = N[i - 1]
                continue
            axis = np.cross(t_prev, t_curr)
            axis_n = np.linalg.norm(axis)
            if axis_n < 1e-9:
                # Antiparallel tangents (curve folds back) — keep N as-is.
                N[i] = N[i - 1]
                continue
            axis = axis / axis_n
            sin_angle = np.sqrt(max(1.0 - cos_angle * cos_angle, 0.0))
            # Rodrigues: v' = v cos θ + (k × v) sin θ + k (k·v) (1 − cos θ)
            v = N[i - 1]
            v_rot = (
                v * cos_angle
                + np.cross(axis, v) * sin_angle
                + axis * np.dot(axis, v) * (1.0 - cos_angle)
            )
            v_rot = v_rot - np.dot(v_rot, t_curr) * t_curr  # re-project
            v_rot = v_rot / np.linalg.norm(v_rot)
            N[i] = v_rot

        # Closed-loop holonomy correction: if the curve closes, distribute
        # the residual twist between N[-1] (transported) and the seed N[0].
        is_closed = bool(np.linalg.norm(pts[-1] - pts[0]) < 1e-6)
        if is_closed:
            holonomy = float(
                np.arctan2(
                    float(np.dot(np.cross(N[-1], N[0]), T[-1])),
                    float(np.dot(N[-1], N[0])),
                ),
            )
            # Vectorised Rodrigues de-twist: each sample's N rotates around
            # its own T by an angle proportional to its arc-length fraction.
            ang = -holonomy * self._curve_s / max(self._curve_total_length, 1e-12)
            cos_a = np.cos(ang)[:, None]
            sin_a = np.sin(ang)[:, None]
            kdotv = np.einsum("ij,ij->i", T, N)[:, None]
            N_rot = N * cos_a + np.cross(T, N) * sin_a + T * kdotv * (1.0 - cos_a)
            N = N_rot / np.linalg.norm(N_rot, axis=1, keepdims=True)
        self._curve_N = N

        self._curve_kdtree = cKDTree(pts)

        # Self-intersection guard.
        if m >= 3:
            d2_dt2 = np.diff(T, axis=0)
            kappa = np.linalg.norm(d2_dt2, axis=1) / np.maximum(np.diff(s_cum), 1e-12)
            kappa_max = float(np.max(kappa)) if kappa.size else 0.0
            if kappa_max * self.radial_max >= 1.0:
                logging.warning(
                    "Sweep: max curvature %.3f × radial_max %.3f ≥ 1 — the "
                    "tube envelope self-intersects in some segments; the "
                    "TPMS may be ill-defined there.",
                    kappa_max,
                    self.radial_max,
                )

    # -- Local-frame field -------------------------------------------------

    def _local_coords(
        self: Sweep,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> tuple[
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
    ]:
        """Return ``(s, r, θ)`` for each input point ``p = (x, y, z)``."""
        shape = x.shape
        pts = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
        _dist, idx = self._curve_kdtree.query(pts, k=1)

        c = self.curve[idx]
        T = self._curve_T[idx]
        N = self._curve_N[idx]
        B = np.cross(T, N)

        delta = pts - c
        s = self._curve_s[idx]
        # signed distance along N and B → radial r and angle θ.
        u = (delta * N).sum(axis=1)
        v = (delta * B).sum(axis=1)
        r = np.sqrt(u * u + v * v)
        theta = np.arctan2(v, u)
        return s.reshape(shape), r.reshape(shape), theta.reshape(shape)

    def _setup_frep_field(self: Sweep) -> None:
        """Build the F-rep field as ``surface_function(k_s s, k_r r, k_θ θ)``."""
        k_x, k_y, k_z = 2.0 * np.pi / self.cell_size
        ps = self.phase_shift

        def _raw_field(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            s, r, theta = self._local_coords(x, y, z)
            return self.surface_function(
                k_x * (s + ps[0]),
                k_y * (r + ps[1]),
                k_z * (theta + ps[2]),
            )

        # Bounds: curve AABB inflated by ``radial_max``.
        bb_min = self.curve.min(axis=0) - self.radial_max
        bb_max = self.curve.max(axis=0) + self.radial_max
        bounds: BoundsType = (
            float(bb_min[0]),
            float(bb_max[0]),
            float(bb_min[1]),
            float(bb_max[1]),
            float(bb_min[2]),
            float(bb_max[2]),
        )

        # Like Conformal, the field is built around discrete data — skip
        # autograd SDF normalisation (it would FD-fall-back and be slow).
        self._raw_field_func = _raw_field
        self._func = _raw_field
        self._bounds = bounds

    # -- Cell-box (tube SDF) -----------------------------------------------

    def _cell_box(self: Sweep) -> Shape:
        """Tube SDF: ``dist_to_curve(p) − radial_max``."""
        from .implicit_ops import from_field

        bounds: BoundsType = (
            self._bounds
            if self._bounds is not None
            else (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
        )

        def _tube_sdf(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            shape = x.shape
            pts = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
            d, _ = self._curve_kdtree.query(pts, k=1)
            return (d - self.radial_max).reshape(shape)

        return from_field(_tube_sdf, bounds=bounds)

    def _isotropic_resolution(self: Sweep) -> int:
        """Use the curve length and tube diameter as the spatial extents."""
        n = int(
            self.resolution
            * max(self._curve_total_length, 2.0 * self.radial_max)
            / max(float(self.cell_size[1]), 1e-12),
        )
        return max(n, 10)

    def _create_grid(
        self: Sweep,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> pv.StructuredGrid:
        """
        Map a parametric (s, r, θ) grid to Cartesian via the curve frames.

        ``x``, ``y``, ``z`` here are the parametric meshgrid coordinates
        from :meth:`Tpms._compute_tpms_field` (they live in (s, r, θ)
        space, *not* Cartesian).  We interpolate the curve point ``c(s)``
        and frame ``(N(s), B(s))`` for each parametric ``s`` value, then
        compute Cartesian ``p = c(s) + r·(cos θ · N + sin θ · B)``.

        The resulting ``StructuredGrid`` has Cartesian point coordinates
        but retains the parametric ordering, so ``self.grid["surface"]``
        (filled by the parent :meth:`Tpms._compute_tpms_field`) can be
        clipped by scalar threshold to extract the TPMS sheet/skeletal —
        same recipe as :class:`CylindricalTpms`.
        """
        # x, y, z come from numpy.meshgrid of parametric linspaces.  Their
        # raw values span [-0.5*cell_size*repeat_cell, +0.5*…] per axis.
        # Map to ``s ∈ [0, total_length]``, ``r ∈ [0, radial_max]``,
        # ``θ ∈ [0, 2π]``.
        s_extent = float(self.cell_size[0]) * float(self.repeat_cell[0])
        s = (x - x.min()) / max(s_extent, 1e-12) * self._curve_total_length
        # Radial: parametric coord shifted to [0, cell_size_r * repeat_cell_r].
        r = y - y.min()
        # Angular: parametric coord scaled to [0, 2π] over the angular
        # extent (cell_size[2] * repeat_cell[2]).
        ang_extent = float(self.cell_size[2]) * float(self.repeat_cell[2])
        theta = (z - z.min()) / max(ang_extent, 1e-12) * 2.0 * np.pi

        # Interpolate curve position + (N, T) frames at each parametric s.
        s_flat = s.ravel(order="F")
        s_clipped = np.clip(s_flat, 0.0, self._curve_total_length)
        c, N, T = (
            _interp_along_curve(s_clipped, self._curve_s, arr)
            for arr in (self.curve, self._curve_N, self._curve_T)
        )
        N /= np.maximum(np.linalg.norm(N, axis=1, keepdims=True), 1e-12)
        T /= np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-12)
        B = np.cross(T, N)

        r_flat = r.ravel(order="F")
        theta_flat = theta.ravel(order="F")
        cos_t = np.cos(theta_flat)
        sin_t = np.sin(theta_flat)

        cart = c + r_flat[:, None] * (cos_t[:, None] * N + sin_t[:, None] * B)
        cart_x = cart[:, 0].reshape(x.shape, order="F")
        cart_y = cart[:, 1].reshape(x.shape, order="F")
        cart_z = cart[:, 2].reshape(x.shape, order="F")

        grid = pv.StructuredGrid(cart_x, cart_y, cart_z)
        grid["coords"] = np.c_[
            x.ravel(order="F"),
            y.ravel(order="F"),
            z.ravel(order="F"),
        ]
        return grid


class Infill(Tpms):
    """Generate a TPMS infill inside a given object."""

    _envelope_mesh_at_full_density = Tpms._envelope_mesh_via_cell_box

    # The parent's ``sheet`` / ``upper_skeletal`` / ``lower_skeletal``
    # properties read from ``self.grid_*`` (legacy grid-clip path), but the
    # density-fitting search optimises against :meth:`generate_vtk` (F-rep
    # marching cubes clipped to the obj envelope) — for small infill objects
    # the two paths discretise to noticeably different volumes.  Re-route
    # these properties to ``generate_vtk`` so ``density=0.5`` reliably gives
    # ``infill.sheet.volume ≈ 0.5 * obj.volume``.
    @property
    def sheet(self: Infill) -> pv.PolyData:
        """Sheet part as a PolyData mesh (uses :meth:`generate_vtk`)."""
        return self.generate_vtk(type_part="sheet")

    @property
    def upper_skeletal(self: Infill) -> pv.PolyData:
        """Upper-skeletal part as a PolyData mesh."""
        return self.generate_vtk(type_part="upper skeletal")

    @property
    def lower_skeletal(self: Infill) -> pv.PolyData:
        """Lower-skeletal part as a PolyData mesh."""
        return self.generate_vtk(type_part="lower skeletal")

    def __init__(
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
        r"""
        Initialize the Infill object.

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
        # Capture the original volume *before* re-orienting normals — for
        # non-manifold inputs (e.g. pyvista's caps-less Cylinder, the
        # Stanford bunny) ``compute_normals(auto_orient_normals=True)`` flips
        # some normals and pyvista's volume sum changes, which would corrupt
        # density-fitting if we relied on ``self.obj.volume`` afterwards.
        self._obj_volume = abs(obj.volume)
        # Auto-orient normals so signed-distance queries (used for the
        # F-rep envelope clip in :meth:`_cell_box`) get the correct sign.
        self.obj = obj.compute_normals(
            auto_orient_normals=True,
            point_normals=True,
            cell_normals=True,
        )
        bounds = np.array(self.obj.bounds)

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
        """
        Density is measured against the input object's volume.

        Uses the volume of the *original* input object rather than the
        re-oriented one stored on ``self.obj`` — see ``__init__``.
        """
        return self._obj_volume

    def _setup_frep_field(self: Infill) -> None:
        """
        Cartesian gyroid field, but bounds set to the obj bbox.

        The plain :class:`Tpms` parent uses origin-centered bounds of size
        ``cell_size * repeat_cell``.  For an Infill of an object that is *not*
        centered at the origin (typical for any imported mesh) those bounds
        miss the object — marching cubes then runs in the wrong place and the
        result isn't actually clipped by the envelope.  Override to use the
        obj bbox so marching cubes covers the right region; the envelope clip
        is enforced in :meth:`_cell_box`.
        """
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

        bounds_arr = np.array(self.obj.bounds, dtype=float)
        bounds: BoundsType = tuple(bounds_arr.tolist())
        self._finalize_frep(_raw_field, bounds)

    def _cell_box(self: Infill) -> Shape:
        """
        Override the cell-box clip with the *object envelope* SDF.

        Without this override the F-rep marching cubes path (used by
        :meth:`Tpms.generate_vtk` and :meth:`Tpms.generate`) would clip the
        TPMS to the cartesian bounding box rather than to the input object —
        producing volumes well above ``obj.volume`` and corrupting density.
        """
        from .implicit_ops import from_field

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


class GradedInfill(Infill):
    """
    Cartesian TPMS infill with offset graded by distance to the envelope.

    Same as :class:`Infill` (TPMS in the world cartesian frame, clipped by
    the input object), except the *thickness* of the TPMS sheet varies with
    the signed distance to the envelope.

    Default behaviour models a **graded shell**: dense material concentrated
    at the skin (``offset_skin = 0.6``, "density 1"), hollowing out toward
    the core (``offset_core = 0.0``, "density 0").  Reverse the parameters
    for a dense-core / porous-skin scaffold.

    The cell pattern is still cartesian (axes-aligned gyroid) — *only* the
    sheet thickness varies in space — so the pattern tiles cleanly even on
    bunny-style envelopes with high-curvature features.
    """

    def __init__(
        self: GradedInfill,
        obj: pv.PolyData,
        surface_function: Field,
        offset_skin: float = 0.6,
        offset_core: float = 0.0,
        transition: float = 0.5,
        smoothness: float = 0.2,
        cell_size: float | Sequence[float] | npt.NDArray[np.float64] | None = None,
        repeat_cell: int | Sequence[int] | npt.NDArray[np.int8] | None = None,
        phase_shift: Sequence[float] = (0.0, 0.0, 0.0),
        resolution: int = 20,
    ) -> None:
        """
        Build a graded TPMS infill.

        :param obj: envelope mesh
        :param surface_function: TPMS function ``f(x,y,z)``
        :param offset_skin: TPMS thickness at the envelope surface
            (``d ≈ 0``).  Default ``0.6`` ⇒ thick / dense skin.
        :param offset_core: TPMS thickness deep in the core (``|d|`` maximal).
            Default ``0.0`` ⇒ no material at the centre (graded shell).
        :param transition: normalised distance ∈ [0, 1] where the offset
            transitions from ``offset_skin`` to ``offset_core``
        :param smoothness: width of the tanh transition (smaller = sharper)
        :param cell_size: unit cell size; mutex with ``repeat_cell``
        :param repeat_cell: number of cells per axis; mutex with ``cell_size``
        :param phase_shift: TPMS phase shift
        :param resolution: per-axis grid resolution
        """
        self._gradation_params = (
            float(offset_skin),
            float(offset_core),
            float(transition),
            float(smoothness),
        )

        graded_offset = self._make_graded_offset_callable(obj, *self._gradation_params)

        super().__init__(
            obj=obj,
            surface_function=surface_function,
            offset=graded_offset,
            cell_size=cell_size,
            repeat_cell=repeat_cell,
            phase_shift=phase_shift,
            resolution=resolution,
            density=None,
        )

    @staticmethod
    def _make_graded_offset_callable(
        envelope: pv.PolyData,
        offset_skin: float,
        offset_core: float,
        transition: float,
        smoothness: float,
    ) -> Field:
        """
        Return ``f(x, y, z) -> offset(p)`` graded by SDF.

        ``d_norm = clip(-d / depth_max, 0, 1)`` ∈ [0, 1] (0 at skin, 1 in
        deepest interior point).  Offset is interpolated from
        ``offset_skin`` (at ``d_norm = 0``) to ``offset_core``
        (at ``d_norm = 1``) via a ``tanh`` profile centred at ``transition``.
        """
        envelope_oriented = envelope.compute_normals(
            auto_orient_normals=True,
            point_normals=True,
            cell_normals=True,
        )

        bounds = np.array(envelope_oriented.bounds, dtype=float)
        coarse = pv.ImageData(
            dimensions=(20, 20, 20),
            spacing=((bounds[1::2] - bounds[::2]) / 19).tolist(),
            origin=bounds[::2].tolist(),
        )
        coarse.compute_implicit_distance(envelope_oriented, inplace=True)
        depth_max = float(max(-np.asarray(coarse["implicit_distance"]).min(), 1e-12))

        def _graded_offset(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            shape = x.shape
            pts = pv.PolyData(np.column_stack([x.ravel(), y.ravel(), z.ravel()]))
            pts.compute_implicit_distance(envelope_oriented, inplace=True)
            d = np.asarray(pts["implicit_distance"])
            d_norm = np.clip(-d / depth_max, 0.0, 1.0)
            # ``sigmoid`` runs 0 → 1 as we go from skin → core.  Multiply that
            # by ``(offset_core - offset_skin)`` and add ``offset_skin`` to
            # interpolate the right way around.
            sigmoid = 0.5 * (
                1.0 + np.tanh((d_norm - transition) / max(smoothness, 1e-6))
            )
            offset_at_p = offset_skin + (offset_core - offset_skin) * sigmoid
            return offset_at_p.reshape(shape)

        return _graded_offset


# Re-export for backward compatibility
from .shape import BoundsType, ShellCreationError  # noqa: E402

__all__ = [
    "CylindricalTpms",
    "GradedInfill",
    "Infill",
    "ShellCreationError",
    "SphericalTpms",
    "Sweep",
    "Tpms",
]
