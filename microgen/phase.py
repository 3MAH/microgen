"""Phase 2.0 — implicit-first, CAD-optional, ``Piece``-aware container.

A :class:`Phase` is a region of space identified by an implicit scalar
field (the canonical representation), with derived materialisations on
demand:

- :meth:`grid` / :meth:`surface_mesh` / :meth:`volume_mesh` — PyVista views
- :attr:`cad` — OCCT BREP via :func:`microgen.cad.shape_to_cad`
- :attr:`pieces` — connected components of ``{field < iso}`` (the
  "Phase = collection of cut/split sub-solids" invariant)
- :attr:`center_of_mass` / :attr:`inertia_matrix` — moments via grid
  quadrature (field-backed) or BRepGProp (CAD-backed)

Three construction paths:

- :class:`Phase` ``(field=..., bounds=..., iso=..., period=...)`` —
  field-first (no CAD required).
- :meth:`Phase.from_shape` — sugar over the field-first path; bridges
  from a :class:`~microgen.shape.shape.Shape`.
- :meth:`Phase.from_cad` — CAD-backed; required when loading a STEP file
  or wrapping a pre-built BREP.

The CAD seam is isolated in :meth:`_materialise_cad`. Switching from
``microgen.cad`` to ``pyvista-cad`` later means swapping that one method;
no other code in microgen needs to change.

Transforms follow the PyVista convention — :meth:`translate`,
:meth:`rotate` and :meth:`scale` take ``inplace=False`` (default; returns
a new :class:`Phase`) or ``inplace=True`` (mutates ``self``); both return
the phase.  :meth:`rotate` and :meth:`scale` pivot on ``point``, the world
origin by default.  :meth:`tile` is a producer (always a new instance).
Every live representation (field, ``iso``, ``period``, CAD, surface mesh,
cached grid) is transformed together so the result stays coherent; derived
caches are invalidated.
"""

from __future__ import annotations

import copy as _copy
import itertools
from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING, Any

import numpy as np

from .shape import _affine

if TYPE_CHECKING:
    from collections.abc import Sequence

    import numpy.typing as npt
    import pyvista as pv
    from scipy.spatial.transform import Rotation

    from .rve import Rve
    from .shape._types import BoundsType, Field, PeriodType
    from .shape.shape import Shape


# Module-level counter for auto-naming (replaces the old mutable
# ``Phase.num_instances`` class attribute, which contaminated test runs).
_PHASE_AUTONAME_COUNTER = itertools.count()


_IMPLICIT_SCALAR = "implicit"


def _bounds_from(lo: npt.ArrayLike, hi: npt.ArrayLike) -> BoundsType:
    """Return ``(xmin, xmax, ymin, ymax, zmin, zmax)`` from two corners."""
    return tuple(float(v) for pair in zip(lo, hi, strict=True) for v in pair)  # type: ignore[return-value]


def _grid_frame(
    sg: pv.StructuredGrid,
) -> tuple[tuple[int, int, int], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Return the dimensions, ``(nx, ny, nz, 3)`` points and edge vectors of a grid.

    The rows of the edge matrix are the steps along ``i``, ``j``, ``k``, so
    node ``(i, j, k)`` sits at ``origin + (i, j, k) @ edges`` for any regular
    grid, including a rotated or scaled one.  A single-node axis gets a unit
    step.
    """
    nx, ny, nz = sg.dimensions
    grid_pts = np.asarray(sg.points, dtype=np.float64).reshape(
        (nx, ny, nz, 3), order="F"
    )
    edges = np.eye(3)
    for axis, n in enumerate((nx, ny, nz)):
        if n > 1:
            index = [0, 0, 0]
            index[axis] = 1
            edges[axis] = grid_pts[tuple(index)] - grid_pts[0, 0, 0]
    return (nx, ny, nz), grid_pts, edges


@dataclass(frozen=True)
class Piece:
    """One connected sub-region of a :class:`Phase`.

    Pieces are what survives "split this phase under periodicity" or
    "raster this phase into a per-cell grid" — they expose the per-piece
    geometric moments without forcing every Phase to be a list of solids.

    Payload fields are populated lazily and may be ``None`` depending on
    how the parent :class:`Phase` was constructed: a CAD-backed phase
    populates ``cad``; a field-backed phase populates ``voxel_mask``;
    a mesh-backed phase populates ``mesh``.
    """

    com: tuple[float, float, float]
    volume: float
    bounds: BoundsType
    cad: Any | None = None
    mesh: pv.PolyData | None = None
    voxel_mask: npt.NDArray[np.bool_] | None = None


class Phase:
    """Microstructure phase (implicit-first, CAD-optional).

    A phase is **one** of the three:

    - field-backed: a callable ``field(x,y,z) -> array`` with negative
      values inside, plus an AABB ``bounds`` and an ``iso`` value (the
      solid is ``{p : field(p) < iso}``).
    - mesh-backed: a triangulated surface ``pv.PolyData``.
    - CAD-backed: a CAD shape (``microgen.cad.CadShape`` today, any
      duck-typed CAD object — including future ``pyvista-cad`` shapes —
      tomorrow).

    The other two materialisations are derived on demand.

    :param field: SDF / level-set ``(x, y, z) -> array``; negative inside.
        Required for field-backed construction.
    :param bounds: axis-aligned bbox ``(xmin, xmax, ymin, ymax, zmin, zmax)``
        spanning the field.  Required if ``field`` is set.
    :param iso: iso-value; the solid is ``{p : field(p) < iso}``.
        Defaults to ``0.0``.
    :param period: ``(Lx, Ly, Lz)`` if the field is intrinsically
        periodic.  Optional.
    :param name: phase name (defaults to auto-generated ``Phase_N``).
    :param resolution: default sampling resolution used by lazy
        :meth:`grid` / :meth:`surface_mesh` / :meth:`pieces`.

    Use :meth:`from_shape`, :meth:`from_cad` or :meth:`from_mesh` to
    build from a :class:`~microgen.shape.shape.Shape`, a pre-built CAD
    object, or a triangulated mesh, respectively.
    """

    def __init__(
        self: Phase,
        *,
        field: Field | None = None,
        bounds: BoundsType | None = None,
        iso: float = 0.0,
        period: PeriodType | None = None,
        name: str | None = None,
        resolution: int = 50,
    ) -> None:
        """Initialize the phase (keyword-only; positional args rejected)."""
        if field is not None and bounds is None:
            err_msg = "bounds must be provided when field is set"
            raise ValueError(err_msg)

        self._field: Field | None = field
        self._bounds: BoundsType | None = bounds
        self._iso: float = float(iso)
        self._period: PeriodType | None = period
        self._resolution: int = int(resolution)
        # Internal caches for non-field-backed payloads. Materialisation
        # rules: CAD-backed phases set ``_cad`` at construction; field-backed
        # phases populate it lazily in :attr:`cad`; mesh-backed phases set
        # ``_surface_mesh`` at construction.
        self._cad: Any | None = None
        self._surface_mesh: pv.PolyData | None = None
        # Set by ``from_grid`` to short-circuit lazy field-sampling in
        # :meth:`grid`; ``None`` for shapes where the field is the source of
        # truth and the grid is regenerated per (bounds, resolution) call.
        self._cached_grid: pv.StructuredGrid | None = None

        self.name: str = (
            name if name is not None else f"Phase_{next(_PHASE_AUTONAME_COUNTER)}"
        )

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_shape(
        cls: type[Phase],
        shape: Shape,
        *,
        bounds: BoundsType | None = None,
        iso: float = 0.0,
        name: str | None = None,
        resolution: int = 50,
    ) -> Phase:
        """Construct a field-backed :class:`Phase` from an implicit :class:`Shape`.

        Inherits ``field``, ``bounds``, and ``period`` from the shape.
        ``bounds`` may be overridden (e.g., to clip a periodic field to a
        sub-region).

        :param shape: source implicit shape; its ``field`` must be set.
        :param bounds: override the shape's bounds; defaults to ``shape.bounds``.
        :param iso: iso-value (default ``0.0``).
        :param name: phase name (auto-generated if omitted).
        :param resolution: default sampling resolution.
        """
        if shape.field is None:
            err_msg = "Cannot build Phase from a Shape without an implicit field"
            raise ValueError(err_msg)
        actual_bounds = bounds if bounds is not None else shape.bounds
        if actual_bounds is None:
            err_msg = (
                "Source Shape has no bounds — pass `bounds=` "
                "explicitly to Phase.from_shape()"
            )
            raise ValueError(err_msg)
        return cls(
            field=shape.field,
            bounds=actual_bounds,
            iso=iso,
            period=shape.period,
            name=name,
            resolution=resolution,
        )

    @classmethod
    def from_cad(
        cls: type[Phase],
        cad: Any,
        *,
        name: str | None = None,
    ) -> Phase:
        """Construct a CAD-backed :class:`Phase`.

        ``cad`` may be a :class:`microgen.cad.CadShape` or any duck-typed
        object with ``.solids()``, ``.center()``, ``.volume()``,
        ``.bounding_box()`` methods (this is the seam that lets a future
        ``pyvista-cad`` backend drop in without touching :class:`Phase`).

        :param cad: a CAD shape (``CadShape`` today)
        :param name: phase name (auto-generated if omitted)
        """
        from .cad import CadShape  # noqa: PLC0415

        if not isinstance(cad, CadShape) and hasattr(cad, "wrapped"):
            cad = CadShape(cad.wrapped)

        instance = cls(name=name)
        instance._cad = cad  # noqa: SLF001
        return instance

    @classmethod
    def from_mesh(
        cls: type[Phase],
        mesh: pv.PolyData,
        *,
        name: str | None = None,
    ) -> Phase:
        """Construct a mesh-backed :class:`Phase` from a closed surface mesh.

        :param mesh: closed triangulated surface
        :param name: phase name (auto-generated if omitted)
        """
        instance = cls(name=name)
        instance._surface_mesh = mesh  # noqa: SLF001
        return instance

    @classmethod
    def from_implicit(
        cls: type[Phase],
        field: Field,
        rve: Rve,
        *,
        iso: float = 0.0,
        period: PeriodType | None = None,
        name: str | None = None,
        resolution: int = 50,
    ) -> Phase:
        """Construct a field-backed :class:`Phase` from a callable + :class:`Rve`.

        Sugar over the field-first constructor that derives ``bounds``
        from the RVE bounding box.

        :param field: implicit scalar field ``(x, y, z) -> array``
            (negative inside).
        :param rve: domain whose AABB becomes the Phase ``bounds``.
        :param iso: iso-value (default ``0.0``).
        :param period: ``(Lx, Ly, Lz)`` if ``field`` is intrinsically periodic.
        :param name: phase name (auto-generated if omitted).
        :param resolution: default sampling resolution.
        """
        bounds = (
            float(rve.min_point[0]),
            float(rve.max_point[0]),
            float(rve.min_point[1]),
            float(rve.max_point[1]),
            float(rve.min_point[2]),
            float(rve.max_point[2]),
        )
        return cls(
            field=field,
            bounds=bounds,
            iso=iso,
            period=period,
            name=name,
            resolution=resolution,
        )

    @classmethod
    def from_grid(
        cls: type[Phase],
        grid: pv.StructuredGrid,
        *,
        scalars: str = "implicit",
        iso: float = 0.0,
        name: str | None = None,
    ) -> Phase:
        """Construct a field-backed :class:`Phase` from a pre-sampled grid.

        Useful when the implicit field is expensive to evaluate (GRF,
        FFT-based fields) and the caller already has a
        :class:`pyvista.StructuredGrid` whose points carry the scalar
        sample. The Phase wraps the grid as its native representation;
        :attr:`grid` returns it untouched, and downstream operations
        (``pieces``, ``surface_mesh``, ``volume_mesh``,
        ``center_of_mass``) read from it directly.

        The Phase's ``field`` is a nearest-neighbour lookup against the
        grid samples — usable for ``Phase.from_shape``-style composition,
        but inexact between sample points.

        :param grid: structured grid whose point data contains the scalar
            field samples.
        :param scalars: name of the scalar array on the grid; copied to
            ``implicit`` on the cached grid when it has another name.
        :param iso: iso-value (default ``0.0``).
        :param name: phase name (auto-generated if omitted).
        """
        if scalars not in grid.point_data:
            err_msg = (
                f"Grid has no point scalar named {scalars!r}. "
                f"Available: {list(grid.point_data.keys())}"
            )
            raise ValueError(err_msg)

        dims, grid_pts, edges = _grid_frame(grid)
        scalar = np.asarray(grid[scalars]).reshape(dims, order="F")
        origin = grid_pts[0, 0, 0].copy()  # a copy: don't keep grid.points alive
        to_index = np.linalg.inv(edges)
        upper = np.array(dims) - 1
        # A regular grid's extreme points are among its 8 corners.
        corners = grid_pts[np.ix_([0, -1], [0, -1], [0, -1])].reshape(-1, 3)
        bounds = _bounds_from(corners.min(axis=0), corners.max(axis=0))

        def _nearest_field(
            x: np.ndarray,
            y: np.ndarray,
            z: np.ndarray,
        ) -> np.ndarray:
            # Node (i, j, k) sits at origin + (i, j, k) @ edges.
            xyz = np.stack(np.broadcast_arrays(x, y, z), axis=-1) - origin
            idx = np.clip(np.round(xyz @ to_index).astype(int), 0, upper)
            return scalar[idx[..., 0], idx[..., 1], idx[..., 2]]

        instance = cls(
            field=_nearest_field,
            bounds=bounds,
            iso=iso,
            name=name,
            resolution=max(dims),
        )
        # Seed the grid cache so subsequent .grid() returns the original
        # (avoids resampling the field via nearest-neighbour).  Every reader
        # (moments, surface mesh, transforms) uses the ``implicit`` array.
        if scalars != _IMPLICIT_SCALAR:
            grid = grid.copy()
            grid[_IMPLICIT_SCALAR] = np.asarray(grid[scalars])
        instance._cached_grid = grid  # noqa: SLF001
        return instance

    # ------------------------------------------------------------------
    # Read-only accessors
    # ------------------------------------------------------------------

    @property
    def field(self: Phase) -> Field | None:
        """The implicit scalar field, or ``None`` for non-field-backed phases."""
        return self._field

    @property
    def bounds(self: Phase) -> BoundsType | None:
        """AABB ``(xmin, xmax, ymin, ymax, zmin, zmax)``.

        For field-backed phases this is set at construction.  For CAD- or
        mesh-backed phases it's derived from the underlying representation.
        """
        if self._bounds is not None:
            return self._bounds
        if self._cad is not None:
            bb = self._cad.bounding_box()
            return (bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax)
        if self._surface_mesh is not None:
            xmin, xmax, ymin, ymax, zmin, zmax = self._surface_mesh.bounds
            return (
                float(xmin),
                float(xmax),
                float(ymin),
                float(ymax),
                float(zmin),
                float(zmax),
            )
        return None

    @property
    def iso(self: Phase) -> float:
        """Iso-value for the implicit field (the solid is ``{field < iso}``)."""
        return self._iso

    @property
    def period(self: Phase) -> PeriodType | None:
        """Intrinsic period ``(Lx, Ly, Lz)`` if the field is periodic."""
        return self._period

    @property
    def resolution(self: Phase) -> int:
        """Default sampling resolution for lazy grid/mesh/pieces materialisation."""
        return self._resolution

    @property
    def is_empty(self: Phase) -> bool:
        """True if the phase has no backing representation."""
        return self._field is None and self._cad is None and self._surface_mesh is None

    # ------------------------------------------------------------------
    # CAD materialisation (the pyvista-cad seam lives here)
    # ------------------------------------------------------------------

    def _materialise_cad(self: Phase) -> Any:
        """Build a CAD representation from the field (or surface mesh).

        This is the **single seam** where :class:`Phase` talks to a CAD
        backend.  Swapping ``microgen.cad`` for ``pyvista-cad`` later
        means rewriting this method only.

        Requires the optional ``[cad]`` extra today.
        """
        if self._field is not None and self._bounds is not None:
            # Wrap field+bounds in a transient Shape and call shape_to_cad.
            from .cad import shape_to_cad  # noqa: PLC0415
            from .shape.shape import Shape  # noqa: PLC0415

            transient = Shape(field=self._field, bounds=self._bounds)
            return shape_to_cad(
                transient, bounds=self._bounds, resolution=self._resolution
            )
        if self._surface_mesh is not None:
            from .cad import mesh_to_shape  # noqa: PLC0415

            mesh = self._surface_mesh
            if not mesh.is_all_triangles:
                mesh = mesh.triangulate()
            triangles = mesh.faces.reshape(-1, 4)[:, 1:]
            points = np.asarray(mesh.points, dtype=np.float64)
            return mesh_to_shape(points, triangles)
        err_msg = "Cannot materialise CAD: phase has no field or mesh"
        raise ValueError(err_msg)

    @cached_property
    def cad(self: Phase) -> Any:
        """The CAD representation (``microgen.cad.CadShape`` today, lazy).

        For CAD-backed phases this returns the stored shape.  For
        field-backed or mesh-backed phases it's lazily materialised via
        :meth:`_materialise_cad`.

        Requires the ``[cad]`` extra (raises ``ImportError`` otherwise).
        """
        if self._cad is not None:
            return self._cad
        return self._materialise_cad()

    # ------------------------------------------------------------------
    # PyVista views (lazy, per-resolution caches via @cached_property are
    # not used here because the resolution kwarg differs per call)
    # ------------------------------------------------------------------

    def grid(self: Phase, resolution: int | None = None) -> pv.StructuredGrid:
        """Return a structured grid sampling of the field.

        Only meaningful for field-backed phases. When the Phase was built
        via :meth:`from_grid`, the cached input grid is returned untouched
        (regardless of the ``resolution`` argument).
        """
        if self._cached_grid is not None:
            return self._cached_grid
        if self._field is None or self._bounds is None:
            err_msg = "grid() requires a field-backed Phase"
            raise ValueError(err_msg)
        import pyvista as pv  # noqa: PLC0415

        res = int(resolution) if resolution is not None else self._resolution
        xmin, xmax, ymin, ymax, zmin, zmax = self._bounds
        xi = np.linspace(xmin, xmax, res)
        yi = np.linspace(ymin, ymax, res)
        zi = np.linspace(zmin, zmax, res)
        x, y, z = np.meshgrid(xi, yi, zi, indexing="ij")
        sg = pv.StructuredGrid(x, y, z)
        sg[_IMPLICIT_SCALAR] = self._field(
            x.ravel(order="F"), y.ravel(order="F"), z.ravel(order="F")
        )
        return sg

    def surface_mesh(self: Phase, resolution: int | None = None) -> pv.PolyData:
        """Return a triangulated surface mesh of the solid boundary.

        - Mesh-backed phase: returns the stored mesh.
        - Field-backed phase: marching cubes on the sampled grid.
        - CAD-backed phase: tessellates via OCCT incremental mesh.
        """
        import pyvista as pv  # noqa: PLC0415

        if self._surface_mesh is not None:
            return self._surface_mesh
        if self._field is not None:
            sg = self.grid(resolution)
            iso = pv.PolyData(
                sg.contour(isosurfaces=[self._iso], scalars=_IMPLICIT_SCALAR)
            )
            return iso.clean().triangulate() if iso.n_cells > 0 else pv.PolyData()
        if self._cad is not None:
            err_msg = (
                "surface_mesh() on a CAD-backed Phase is not implemented yet "
                "(would need OCCT BRepMesh_IncrementalMesh tessellation)."
            )
            raise NotImplementedError(err_msg)
        err_msg = "Cannot build surface_mesh on an empty Phase"
        raise ValueError(err_msg)

    def volume_mesh(self: Phase, resolution: int | None = None) -> pv.UnstructuredGrid:
        """Return the volumetric cells where ``field < iso``.

        Only meaningful for field-backed phases.
        """
        if self._field is None:
            err_msg = "volume_mesh() requires a field-backed Phase"
            raise ValueError(err_msg)
        sg = self.grid(resolution)
        return sg.clip_scalar(scalars=_IMPLICIT_SCALAR, value=self._iso, invert=True)

    # ------------------------------------------------------------------
    # Pieces — the "Phase = collection of cut/split sub-solids" invariant
    # ------------------------------------------------------------------

    @cached_property
    def pieces(self: Phase) -> list[Piece]:
        """Connected components of the phase (the "sub-pieces" invariant).

        - CAD-backed phase: one :class:`Piece` per ``TopoDS_Solid``.
        - Field-backed phase: ``scipy.ndimage.label`` on
          ``grid_field < iso``; one piece per connected component.
        - Mesh-backed phase: ``polydata.connectivity().split_bodies()``.
        """
        if self._cad is not None:
            return self._pieces_from_cad()
        if self._field is not None:
            return self._pieces_from_field()
        if self._surface_mesh is not None:
            return self._pieces_from_mesh()
        return []

    def _pieces_from_cad(self: Phase) -> list[Piece]:
        from .cad import CadShape  # noqa: PLC0415

        out: list[Piece] = []
        for solid in self._cad.solids():
            wrapped = solid if isinstance(solid, CadShape) else CadShape(solid)
            c = wrapped.center()
            bb = wrapped.bounding_box()
            out.append(
                Piece(
                    com=(float(c.x), float(c.y), float(c.z)),
                    volume=float(wrapped.volume()),
                    bounds=(bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax),
                    cad=wrapped,
                )
            )
        return out

    def _pieces_from_field(self: Phase) -> list[Piece]:
        from scipy.ndimage import (  # noqa: PLC0415
            center_of_mass,
            find_objects,
            label,
        )

        sg = self.grid()
        dims, grid_pts, edges = _grid_frame(sg)
        scalar = np.asarray(sg[_IMPLICIT_SCALAR]).reshape(dims, order="F")
        inside = scalar < self._iso
        labels, n_labels = label(inside)
        if n_labels == 0:
            return []

        cell_volume = float(abs(np.linalg.det(edges)))
        origin = grid_pts[0, 0, 0]
        coms = center_of_mass(inside, labels=labels, index=range(1, n_labels + 1))
        slices = find_objects(labels)

        out: list[Piece] = []
        for label_id in range(1, n_labels + 1):
            sl = slices[label_id - 1]
            local = labels[sl] == label_id
            mask = np.zeros(dims, dtype=bool)
            mask[sl] = local
            voxel_count = int(local.sum())
            com = origin + np.asarray(coms[label_id - 1]) @ edges
            com_world = (float(com[0]), float(com[1]), float(com[2]))
            nodes = grid_pts[sl][local]
            piece_bounds = _bounds_from(nodes.min(axis=0), nodes.max(axis=0))
            out.append(
                Piece(
                    com=com_world,
                    volume=voxel_count * cell_volume,
                    bounds=piece_bounds,
                    voxel_mask=mask,
                )
            )
        return out

    def _pieces_from_mesh(self: Phase) -> list[Piece]:
        out: list[Piece] = []
        for body in self._surface_mesh.split_bodies():  # type: ignore[union-attr]
            poly = body.extract_surface()
            com = poly.center_of_mass()
            xmin, xmax, ymin, ymax, zmin, zmax = poly.bounds
            out.append(
                Piece(
                    com=(float(com[0]), float(com[1]), float(com[2])),
                    volume=float(poly.volume),
                    bounds=(
                        float(xmin),
                        float(xmax),
                        float(ymin),
                        float(ymax),
                        float(zmin),
                        float(zmax),
                    ),
                    mesh=poly,
                )
            )
        return out

    # ------------------------------------------------------------------
    # Moments
    # ------------------------------------------------------------------

    @cached_property
    def center_of_mass(self: Phase) -> npt.NDArray[np.float64]:
        """Volumetric center of mass.

        For field-backed phases this is computed by quadrature on the
        sampled grid (no OCCT needed).  For CAD-backed phases this uses
        ``BRepGProp.VolumeProperties_s``.
        """
        if self._cad is not None:
            from OCP.BRepGProp import BRepGProp  # noqa: PLC0415
            from OCP.GProp import GProp_GProps  # noqa: PLC0415

            props = GProp_GProps()
            BRepGProp.VolumeProperties_s(self._cad.wrapped, props)
            com = props.CentreOfMass()
            return np.array([com.X(), com.Y(), com.Z()])
        if self._field is not None:
            pts, _ = self._inside_nodes("center_of_mass")
            return pts.mean(axis=0)
        err_msg = "Cannot compute center_of_mass on an empty Phase"
        raise ValueError(err_msg)

    def _inside_nodes(
        self: Phase, caller: str
    ) -> tuple[npt.NDArray[np.float64], float]:
        """Return the grid nodes inside the solid and the volume of one grid cell.

        The cell volume is ``|det(e_i, e_j, e_k)|`` of the grid's edge vectors,
        so it stays right for a non-cubic, rotated or scaled cached grid.
        """
        sg = self.grid()
        inside = np.asarray(sg[_IMPLICIT_SCALAR]) < self._iso
        if not inside.any():
            err_msg = f"{caller}: field is positive everywhere on the grid"
            raise ValueError(err_msg)
        _, _, edges = _grid_frame(sg)
        pts = np.asarray(sg.points, dtype=np.float64)
        return pts[inside], float(abs(np.linalg.det(edges)))

    @cached_property
    def inertia_matrix(self: Phase) -> npt.NDArray[np.float64]:
        """Inertia tensor of the phase about its center of mass (unit density).

        Field-backed: grid quadrature.  CAD-backed: ``BRepGProp``.  Both
        use the center of mass as reference point, so the tensor is
        invariant under :meth:`translate`.
        """
        if self._cad is not None:
            from OCP.BRepGProp import BRepGProp  # noqa: PLC0415
            from OCP.GProp import GProp_GProps  # noqa: PLC0415

            props = GProp_GProps()
            BRepGProp.VolumeProperties_s(self._cad.wrapped, props)
            inm = props.MatrixOfInertia()
            return np.array(
                [
                    [inm.Value(1, 1), inm.Value(1, 2), inm.Value(1, 3)],
                    [inm.Value(2, 1), inm.Value(2, 2), inm.Value(2, 3)],
                    [inm.Value(3, 1), inm.Value(3, 2), inm.Value(3, 3)],
                ]
            )
        if self._field is not None:
            pts, cell_volume = self._inside_nodes("inertia_matrix")
            x, y, z = (pts - self.center_of_mass).T
            ixx = float(((y * y + z * z) * cell_volume).sum())
            iyy = float(((x * x + z * z) * cell_volume).sum())
            izz = float(((x * x + y * y) * cell_volume).sum())
            ixy = -float((x * y * cell_volume).sum())
            ixz = -float((x * z * cell_volume).sum())
            iyz = -float((y * z * cell_volume).sum())
            return np.array(
                [[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]],
                dtype=np.float64,
            )
        err_msg = "Cannot compute inertia_matrix on an empty Phase"
        raise ValueError(err_msg)

    # ------------------------------------------------------------------
    # Copy
    # ------------------------------------------------------------------

    def copy(self: Phase) -> Phase:
        """Return an independent copy with every live representation duplicated.

        The field callable and (immutable) bounds tuple are shared; the CAD
        (:meth:`CadShape.copy`, i.e. ``BRepBuilderAPI_Copy``), the surface mesh
        and any cached grid are deep-copied, so an ``inplace=False``
        transform on the copy never touches the original.
        """
        new = self._shallow_copy()
        if self._cad is not None:
            new._cad = self._cad.copy()  # noqa: SLF001
        if self._surface_mesh is not None:
            new._surface_mesh = self._surface_mesh.copy()  # noqa: SLF001
        if self._cached_grid is not None:
            new._cached_grid = self._cached_grid.copy()  # noqa: SLF001
        return new

    def _invalidate_derived(self: Phase) -> None:
        """Drop cached ``@cached_property`` results after a mutation."""
        for key in ("cad", "pieces", "center_of_mass", "inertia_matrix"):
            self.__dict__.pop(key, None)

    # ------------------------------------------------------------------
    # Transforms — PyVista convention: imperative verb + ``inplace=False``.
    #
    # ``inplace=False`` (default) returns a new Phase; ``inplace=True``
    # mutates ``self``.  Both return the Phase.  ``_apply_affine`` maps
    # every live representation (field, iso, period, CAD, surface mesh,
    # cached grid) through the same affine map, so the result stays
    # coherent whichever mode the caller picked.
    # ------------------------------------------------------------------

    def translate(
        self: Phase, offset: Sequence[float], *, inplace: bool = False
    ) -> Phase:
        """Translate the phase by ``offset``.

        ``period`` and ``iso`` are kept.

        :param offset: ``(dx, dy, dz)`` world-space shift.
        :param inplace: mutate ``self`` when ``True``; otherwise (default)
            return a new translated :class:`Phase`.
        :return: the translated phase.
        :raises ValueError: if the phase is empty.
        """
        return self._transform(_affine.translation(offset), "translate", inplace)

    def rotate(
        self: Phase,
        rotation: Rotation | npt.ArrayLike,
        point: Sequence[float] | None = None,
        *,
        inplace: bool = False,
    ) -> Phase:
        """Rotate the phase about ``point``.

        Field, CAD, surface mesh and cached grid rotate together.  A
        rotation within ``1e-12`` of a signed permutation (quarter turns
        about the axes) is snapped to it and permutes ``period``; any other
        rotation sets ``period`` to ``None`` and clips the field to the
        rotated original ``bounds``, so the solid sampled on the new,
        larger AABB is exactly the rotated solid (no extra material from a
        periodic or unbounded field).

        :param rotation: a SciPy ``Rotation`` or a proper orthogonal 3x3
            matrix, following PyVista's ``RotationLike``.
        :param point: pivot; defaults to the world origin.
        :param inplace: mutate ``self`` when ``True``; otherwise (default)
            return a new rotated :class:`Phase`.
        :return: the rotated phase.
        :raises ValueError: if the phase is empty or the matrix is not a
            rotation.
        """
        return self._transform(
            _affine.rotation_about(rotation, point), "rotate", inplace
        )

    def scale(
        self: Phase,
        factor: float | Sequence[float],
        point: Sequence[float] | None = None,
        *,
        inplace: bool = False,
    ) -> Phase:
        r"""Scale the phase about ``point``.

        The field becomes :math:`c\, f(p + (x - p) / s)` with
        :math:`c = \min_i s_i` and ``iso`` is multiplied by :math:`c`, so the
        solid is exactly the scaled solid; a signed distance stays exact for
        a uniform factor and becomes a lower bound of the distance for a
        per-axis one.  ``period`` scales with the factors.  Pass
        ``point=phase.center_of_mass`` to scale about the phase itself.

        :param factor: uniform factor or ``(sx, sy, sz)``, strictly positive.
        :param point: pivot; defaults to the world origin.
        :param inplace: mutate ``self`` when ``True``; otherwise (default)
            return a new scaled :class:`Phase`.
        :return: the scaled phase.
        :raises ValueError: if the phase is empty or a factor is not
            strictly positive.
        """
        return self._transform(_affine.scaling(factor, point), "scale", inplace)

    def _transform(
        self: Phase, affine: _affine.AffineMap, verb: str, inplace: bool
    ) -> Phase:
        if self.is_empty:
            err_msg = f"Cannot {verb} an empty Phase"
            raise ValueError(err_msg)
        target = self if inplace else self._shallow_copy()
        target._apply_affine(affine)  # noqa: SLF001
        return target

    def _shallow_copy(self: Phase) -> Phase:
        """Copy sharing every representation; ``_apply_affine`` replaces them all."""
        new = _copy.copy(self)
        new._invalidate_derived()  # noqa: SLF001
        return new

    def _apply_affine(self: Phase, affine: _affine.AffineMap) -> None:
        """Map every live representation through ``affine`` (new objects, no mutation)."""
        if self._field is not None and self._bounds is not None:
            field = self._field
            if affine.grows_bounds:
                field = self._clipped_to_bounds(field)
            self._field = affine.apply_field(field)
            self._bounds = affine.apply_bounds(self._bounds)
            self._iso *= affine.multiplier
        self._period = affine.apply_period(self._period)
        if self._cad is not None:
            from .cad import transform_geometry  # noqa: PLC0415

            self._cad = transform_geometry(self._cad, affine.homogeneous[:3])
        if self._surface_mesh is not None:
            self._surface_mesh = self._surface_mesh.transform(
                affine.homogeneous, inplace=False
            )
        if self._cached_grid is not None:
            grid = self._cached_grid.transform(affine.homogeneous, inplace=False)
            if affine.multiplier != 1.0:
                grid[_IMPLICIT_SCALAR] = affine.multiplier * np.asarray(
                    grid[_IMPLICIT_SCALAR]
                )
            self._cached_grid = grid
        self._invalidate_derived()

    def _clipped_to_bounds(self: Phase, field: Field) -> Field:
        r"""Restrict the solid ``{field < iso}`` to the current ``bounds`` box.

        Returns :math:`\max(f, d_{box} + \mathrm{iso})`, whose
        ``iso``-level set is the solid intersected with the box.
        """
        from .shape.implicit_ops import box  # noqa: PLC0415

        b = np.asarray(self._bounds, dtype=np.float64).reshape(3, 2)
        box_field = box(dims=b[:, 1] - b[:, 0], center=b.mean(axis=1)).require_field()
        iso = self._iso

        def _clipped(x, y, z):  # noqa: ANN001, ANN202
            return np.maximum(field(x, y, z), box_field(x, y, z) + iso)

        return _clipped

    def tile(self: Phase, rve: Rve, grid: tuple[int, int, int]) -> Phase:
        r"""Periodically tile the phase on the RVE over ``grid`` copies.

        A **producer**: always returns a new :class:`Phase` (no ``inplace``).
        Builds ``grid[0] * grid[1] * grid[2]`` copies of the current phase,
        translated by :math:`-\mathrm{dim}\,(n/2 - 1/2 - i)` along each axis
        (:math:`i = 0 \dots n-1`), i.e. a supercell centred on the original:
        ``grid=(1, 1, 1)`` leaves it in place, and an even count shifts the
        copies by half a cell.  The copies are fused.  Only implemented for
        CAD-backed phases today (the field-backed equivalent, domain folding
        via ``mod``, lands in a follow-up; the periodic-shape work in
        :mod:`microgen.shape.implicit_ops` already covers it for
        :class:`Shape` directly).

        :param rve: RVE whose ``dim`` sets the tiling step.
        :param grid: number of copies along x, y and z.
        :return: a new CAD-backed :class:`Phase` holding the fused copies.
        :raises NotImplementedError: if the phase is not CAD-backed.
        """
        if self._cad is None:
            err_msg = (
                "Phase.tile is implemented for CAD-backed phases today; "
                "field-backed tiling lives on the source Shape (use "
                "microgen.shape.implicit_ops.repeat there)."
            )
            raise NotImplementedError(err_msg)
        from .cad import (  # noqa: PLC0415
            make_compound_from_solids,
            translate_solid,
        )

        copies = []
        for idx in np.ndindex(*grid):
            # Copies are centred on the original: offsets are multiples of
            # ``rve.dim`` symmetric about zero (``grid=(1, 1, 1)`` is a copy).
            offset = -rve.dim * (0.5 * np.array(grid) - 0.5 - np.array(idx))
            copies.append(translate_solid(self._cad.wrapped, offset))
        new = Phase(name=self.name, resolution=self._resolution)
        new._cad = make_compound_from_solids(copies)  # noqa: SLF001
        return new

    # ------------------------------------------------------------------
    # Misc
    # ------------------------------------------------------------------

    def __repr__(self: Phase) -> str:
        kind = (
            "field"
            if self._field is not None
            else "cad"
            if self._cad is not None
            else "mesh"
            if self._surface_mesh is not None
            else "empty"
        )
        return f"Phase(name={self.name!r}, kind={kind!r}, bounds={self.bounds})"


__all__ = ["Phase", "Piece"]
