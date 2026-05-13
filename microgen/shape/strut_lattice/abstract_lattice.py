"""
=======================================================================
Abstract Lattice (:mod:`microgen.shape.strut_lattice.abstract_lattice`)
=======================================================================
"""

from __future__ import annotations

import sys
import warnings
from abc import abstractmethod
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.optimize import root_scalar
from scipy.spatial.transform import Rotation

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

from ...cad import CadShape
from ...mesh import mesh, mesh_periodic
from ...operations import fuse_shapes
from ...periodic import periodic_split_and_translate
from ...phase import Phase
from ...rve import Rve
from ..box import Box
from ..cylinder import Cylinder
from ..shape import Shape
from ..sphere import Sphere

if TYPE_CHECKING:
    from microgen.shape import KwargsGenerateType, Vector3DType

BALL_POINT_RADIUS_TOLERANCE = 1e-5
_DENSITY_ROOT_RADIUS_MIN = 1e-3
_DEFAULT_STRUT_RADIUS = 0.05


class AbstractLattice(Shape):
    """Abstract class for strut-based lattices, configured via method chaining.

    Every parameter has a default: ``__init__`` is callable with no arguments
    and produces a lattice with ``strut_radius=0.05`` and ``cell_size=1.0``.
    All parameters (radius / density, cell size, strut heights, joints,
    center, orientation, custom vertices / pairs) are also exposed via chained
    ``with_*`` setters; heavy work (vertex layout, density root-finding, CAD
    assembly) runs lazily on the first call to :meth:`generate_cad` /
    :meth:`generate_surface_mesh`.

    Example::

        lattice = OctetTruss().with_cell_size(1.0).generate_surface_mesh()

    Mutual exclusivity: ``strut_radius`` and ``density`` cannot both be set
    explicitly.  Passing both to ``__init__`` raises; the ``with_*`` setters
    emit a :class:`UserWarning` and clear the previous value on a mode switch.
    """

    _UNIT_CUBE_SIZE = 1.0
    _DEFAULT_STRUT_HEIGHTS: float | list[float] | None = None

    def __init__(
        self,
        strut_radius: float | None = None,
        *,
        density: float | None = None,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType | Rotation = (0, 0, 0),
    ) -> None:
        """Initialize with ``strut_radius`` (or ``density``) plus defaults.

        Exactly one of ``strut_radius`` / ``density`` should be set, either
        here or later via the chained ``with_strut_radius`` / ``with_density``
        setters.  The remaining geometry parameters (``cell_size``,
        ``strut_heights``, ``strut_joints``, custom vertices / pairs) are
        configured via the ``with_*`` setters.
        """
        super().__init__(center=center, orientation=orientation)

        if strut_radius is not None and density is not None:
            err_msg = (
                "strut_radius and density cannot both be set. "
                "Pass one positionally and leave the other unset, or use "
                ".with_strut_radius() / .with_density() to choose later."
            )
            raise ValueError(err_msg)
        if density is not None and not 0.0 < density <= 1.0:
            err_msg = f"density must be between 0 and 1. Given: {density}"
            raise ValueError(err_msg)

        # Configuration (mutated by chained setters).
        self._strut_radius: float | None = strut_radius
        self._strut_heights: float | list[float] | None = type(
            self
        )._DEFAULT_STRUT_HEIGHTS
        self._user_base_vertices: npt.NDArray[np.float64] | None = None
        self._user_strut_vertex_pairs: npt.NDArray[np.int64] | None = None
        self._cell_size: float = 1.0
        self._strut_joints: bool = False
        self._density: float | None = density

        # Track which of strut_radius / density was last set by the user, so
        # the setters can warn on a mode switch without false-firing after the
        # lazy density-fit writes back to ``_strut_radius``.
        self._strut_radius_explicit: bool = strut_radius is not None
        self._density_explicit: bool = density is not None

        # Lazy caches (invalidated by setters).
        self._cad_shape: CadShape | None = None
        self._vtk_shape: tuple[tuple[float, int, bool], pv.PolyData] | None = None
        self._geometry: dict | None = None
        self._rve: Rve | None = None

    # ------------------------------------------------------------------
    # Chained setters
    # ------------------------------------------------------------------

    def with_strut_radius(self: Self, radius: float) -> Self:
        """Set the strut radius.

        Clears any previously set density (last-set wins between strut radius
        and density).  Emits a :class:`UserWarning` if density was already
        explicitly set, so the override is visible to the caller.
        """
        if self._density_explicit:
            warnings.warn(
                "Overriding explicit density with strut_radius; the previous "
                "density value will be cleared.",
                stacklevel=2,
            )
            self._density_explicit = False
        self._strut_radius_explicit = True
        self._strut_radius = radius
        self._density = None
        self._invalidate_caches()
        return self

    def with_strut_heights(
        self: Self,
        heights: float | list[float],
    ) -> Self:
        """Set strut heights, as a scalar or per-strut list (unit-cell units)."""
        self._strut_heights = heights
        self._invalidate_caches()
        return self

    def with_cell_size(self: Self, size: float) -> Self:
        """Set the cubic cell edge length."""
        self._cell_size = size
        self._invalidate_caches()
        return self

    def with_strut_joints(self: Self, *, enabled: bool = True) -> Self:
        """Enable (default) or disable spherical joints at vertices."""
        self._strut_joints = enabled
        self._invalidate_caches()
        return self

    def with_density(self: Self, density: float) -> Self:
        """Set target density in (0, 1].

        Clears any previously set strut radius (last-set wins between strut
        radius and density).  Emits a :class:`UserWarning` if strut_radius
        was already explicitly set, so the override is visible to the caller.
        """
        if not 0.0 < density <= 1.0:
            err_msg = f"density must be between 0 and 1. Given: {density}"
            raise ValueError(err_msg)
        if self._strut_radius_explicit:
            warnings.warn(
                "Overriding explicit strut_radius with density; the previous "
                "strut_radius value will be cleared.",
                stacklevel=2,
            )
            self._strut_radius_explicit = False
        self._density_explicit = True
        self._density = density
        self._strut_radius = None
        self._invalidate_caches()
        return self

    def with_base_vertices(
        self: Self,
        vertices: npt.NDArray[np.float64],
    ) -> Self:
        """Override the subclass's default vertex layout (unit-cube units)."""
        self._user_base_vertices = vertices
        self._invalidate_caches()
        return self

    def with_strut_vertex_pairs(
        self: Self,
        pairs: npt.NDArray[np.int64],
    ) -> Self:
        """Override the subclass's default strut connectivity."""
        self._user_strut_vertex_pairs = pairs
        self._invalidate_caches()
        return self

    def with_center(self: Self, center: Vector3DType) -> Self:
        """Set the lattice center; clears caches so the next generate reuses it."""
        self._center = center
        self._invalidate_caches()
        return self

    def with_orientation(self: Self, orientation: Vector3DType | Rotation) -> Self:
        """Set the lattice orientation (Euler ZXZ degrees or :class:`Rotation`)."""
        self._orientation = (
            orientation
            if isinstance(orientation, Rotation)
            else Rotation.from_euler("ZXZ", orientation, degrees=True)
        )
        self._invalidate_caches()
        return self

    def _invalidate_caches(self) -> None:
        self._cad_shape = None
        self._vtk_shape = None
        self._geometry = None
        self._rve = None

    # ------------------------------------------------------------------
    # Public read accessors
    # ------------------------------------------------------------------

    @property
    def cell_size(self) -> float:
        """Cubic cell edge length."""
        return self._cell_size

    @property
    def strut_joints(self) -> bool:
        """Whether spherical joints are added at vertices."""
        return self._strut_joints

    @property
    def density(self) -> float | None:
        """User-requested target density (or ``None`` if radius-driven)."""
        return self._density

    @property
    def strut_radius(self) -> float:
        """Effective strut radius.

        If only :meth:`with_density` was set, this triggers the lazy
        density-to-radius root-find on first access.  If neither
        ``strut_radius`` nor ``density`` was explicitly chosen, falls back to
        :data:`_DEFAULT_STRUT_RADIUS`.
        """
        if self._strut_radius is None and self._density is not None:
            self._strut_radius = self._compute_radius_to_fit_density()
        if self._strut_radius is None:
            self._strut_radius = _DEFAULT_STRUT_RADIUS
        return self._strut_radius

    @property
    def base_vertices(self) -> npt.NDArray[np.float64]:
        """Vertex coordinates in a unit cubic RVE centered on the origin."""
        if self._user_base_vertices is not None:
            return self._user_base_vertices
        return self._generate_base_vertices()

    @property
    def strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        """Pairs of vertex indices defining each strut."""
        if self._user_strut_vertex_pairs is not None:
            return self._user_strut_vertex_pairs
        return self._generate_strut_vertex_pairs()

    @property
    def strut_number(self) -> int:
        """Number of struts in the lattice."""
        return len(self.strut_vertex_pairs)

    @property
    def strut_heights(self) -> list[float]:
        """Per-strut height list, scaled by :attr:`cell_size`."""
        if self._strut_heights is None:
            err_msg = (
                "strut_heights must be defined by the subclass or set "
                "via .with_strut_heights()"
            )
            raise NotImplementedError(err_msg)
        if isinstance(self._strut_heights, float):
            return [self._strut_heights * self._cell_size] * self.strut_number
        return [h * self._cell_size for h in self._strut_heights]

    @property
    def rve(self) -> Rve:
        """RVE that bounds the lattice cell."""
        if self._rve is None:
            self._rve = Rve(dim=self._cell_size, center=self.center)
        return self._rve

    @property
    def vertices(self) -> npt.NDArray[np.float64]:
        """Lattice vertex coordinates in world space."""
        return self._geom()["vertices"]

    @property
    def strut_centers(self) -> npt.NDArray[np.float64]:
        """Midpoints of each strut."""
        return self._geom()["strut_centers"]

    @property
    def strut_directions_cartesian(self) -> npt.NDArray[np.float64]:
        """Unit direction vectors of each strut."""
        return self._geom()["strut_directions"]

    @property
    def strut_rotations(self) -> list[Rotation]:
        """Rotations bringing the default x-axis cylinder onto each strut."""
        return self._geom()["strut_rotations"]

    # ------------------------------------------------------------------
    # Subclass hooks
    # ------------------------------------------------------------------

    @abstractmethod
    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        """Return vertex coordinates in a unit cube centered on the origin."""

    @abstractmethod
    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        """Return pairs of vertex indices defining each strut."""

    # ------------------------------------------------------------------
    # Lazy geometry / validation
    # ------------------------------------------------------------------

    def _geom(self) -> dict:
        if self._geometry is not None:
            return self._geometry
        self._validate()
        vertices = self.center + self._cell_size * self.base_vertices
        pairs = self.strut_vertex_pairs
        strut_centers = np.mean(vertices[pairs], axis=1)
        vectors = np.diff(vertices[pairs], axis=1).squeeze()
        directions = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
        rotations = self._compute_rotations_from(directions)
        self._geometry = {
            "vertices": vertices,
            "strut_centers": strut_centers,
            "strut_directions": directions,
            "strut_rotations": rotations,
        }
        return self._geometry

    def _validate(self) -> None:
        if self._strut_heights is None:
            err_msg = "strut_heights must be defined by the subclass"
            raise NotImplementedError(err_msg)
        if (
            isinstance(self._strut_heights, list)
            and len(self._strut_heights) != self.strut_number
        ):
            err_msg = (
                f"strut_heights must contain {self.strut_number} values, "
                f"but {len(self._strut_heights)} were provided."
            )
            raise ValueError(err_msg)

    @staticmethod
    def _compute_rotations_from(
        directions: npt.NDArray[np.float64],
    ) -> list[Rotation]:
        default_direction = np.array([1.0, 0.0, 0.0])
        rotations_list: list[Rotation] = []
        for i in range(directions.shape[0]):
            d = directions[i]
            if np.all(d == default_direction) or np.all(d == -default_direction):
                rotations_list.append(Rotation.from_rotvec(np.zeros(3)))
            else:
                rotation, _ = Rotation.align_vectors(d, default_direction)
                rotations_list.append(rotation)
        return rotations_list

    def _compute_radius_to_fit_density(self) -> float:
        """Root-find the strut radius that matches :attr:`density`.

        Does NOT cache the CAD produced during the search: ``root_scalar``'s
        final evaluation is at the bracket midpoint, not exactly at the root,
        so the last CAD does not correspond to the returned radius.  The
        caller (``generate_cad``) rebuilds at ``result.root``.
        """

        def calc_density(radius: float) -> float:
            self._strut_radius = radius
            return self._generate_cad().volume() / (self._cell_size**3)

        result = root_scalar(
            lambda r: calc_density(r) - self._density,
            bracket=[_DENSITY_ROOT_RADIUS_MIN, self._cell_size],
        )
        return float(result.root)

    # ------------------------------------------------------------------
    # Terminal generation
    # ------------------------------------------------------------------

    def generate_cad(self, **_: KwargsGenerateType) -> CadShape:
        """Generate (or return cached) CAD shape."""
        if isinstance(self._cad_shape, CadShape):
            return self._cad_shape
        # Trigger lazy radius resolution + validation via the property.
        _ = self.strut_radius
        self._cad_shape = self._generate_cad()
        return self._cad_shape

    cad_shape = property(generate_cad)

    @property
    def volume(self) -> float:
        """Volume of the generated CAD shape."""
        return self.cad_shape.volume()

    def _generate_cad(self, **_: KwargsGenerateType) -> CadShape:
        """Build the CAD lattice from the current configuration."""
        list_phases: list[Phase] = []
        list_periodic_phases: list[Phase] = []

        for i in range(self.strut_number):
            strut = Cylinder(
                center=tuple(self.strut_centers[i]),
                orientation=self.strut_rotations[i],
                height=self.strut_heights[i],
                radius=self._strut_radius,
            )
            list_phases.append(Phase(strut.generate_cad()))
        if self._strut_joints:
            for vertex in self.vertices:
                joint = Sphere(center=tuple(vertex), radius=self._strut_radius)
                list_phases.append(Phase(joint.generate_cad()))

        for phase in list_phases:
            list_periodic_phases.append(
                periodic_split_and_translate(phase=phase, rve=self.rve),
            )

        lattice = fuse_shapes(
            [phase.shape for phase in list_periodic_phases],
            retain_edges=False,
        )

        bounding_box = Box(
            center=self.center,
            dim=(self._cell_size, self._cell_size, self._cell_size),
        ).generate_cad()

        return bounding_box.intersect(lattice)

    def generate_surface_mesh(
        self,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Return a surface mesh of the lattice (for visualisation).

        Delegates to :meth:`mesh_for_fem` with default parameters
        (``size=0.02, order=1, periodic=True``).  Callers who need to
        control mesh size, element order, or periodicity should call
        :meth:`mesh_for_fem` directly.
        """
        return self.mesh_for_fem()

    vtk_shape = property(generate_surface_mesh)

    def mesh_for_fem(
        self,
        size: float = 0.02,
        order: int = 1,
        *,
        periodic: bool = True,
    ) -> pv.PolyData:
        """Generate (or return cached) the lattice FEM tet-mesh surface.

        Path: ``cad_shape`` -> STEP -> gmsh (:func:`microgen.mesh_periodic`
        or :func:`microgen.mesh`) -> ``pv.read`` -> :meth:`extract_surface`.
        Requires the ``[cad]`` extra and gmsh.  Cached per
        ``(size, order, periodic)`` tuple on the instance.
        """
        params = (size, order, periodic)
        if self._vtk_shape is not None:
            cached_params, cached_mesh = self._vtk_shape
            if cached_params == params:
                return cached_mesh

        cad_lattice = self.cad_shape
        list_phases = [Phase(cad_lattice)]

        with (
            NamedTemporaryFile(suffix=".step", delete=False) as cad_step_file,
            NamedTemporaryFile(suffix=".vtk", delete=False) as mesh_file,
        ):
            cad_lattice.export_step(cad_step_file.name)
            mesher = mesh_periodic if periodic else mesh
            mesher_kwargs = (
                {"rve": self.rve, "list_phases": list_phases}
                if periodic
                else {"list_phases": list_phases}
            )
            mesher(
                mesh_file=cad_step_file.name,
                size=size,
                order=order,
                output_file=mesh_file.name,
                **mesher_kwargs,
            )
            vtk_lattice = pv.read(mesh_file.name).extract_surface(algorithm=None)

        # Solve compatibility issues of NamedTemporaryFiles with Windows.
        for file in (cad_step_file.name, mesh_file.name):
            Path(file).unlink()

        self._vtk_shape = (params, vtk_lattice)
        return vtk_lattice
