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
a new :class:`Phase`) or ``inplace=True`` (mutates ``self``).  :meth:`tile`
is a producer (always a new instance).  Every live representation (field,
CAD, surface mesh, cached grid) is transformed together so the result
stays coherent; derived caches are invalidated.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy.spatial.transform import Rotation

if TYPE_CHECKING:
    from collections.abc import Sequence

    import numpy.typing as npt
    import pyvista as pv

    from .rve import Rve
    from .shape._types import BoundsType, Field, PeriodType
    from .shape.shape import Shape


# Module-level counter for auto-naming (replaces the old mutable
# ``Phase.num_instances`` class attribute, which contaminated test runs).
_PHASE_AUTONAME_COUNTER = itertools.count()


_IMPLICIT_SCALAR = "implicit"


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

        :param shape: source implicit shape; must have ``func is not None``.
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
        func: Field,
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

        :param func: implicit scalar field ``(x, y, z) -> array``
            (negative inside).
        :param rve: domain whose AABB becomes the Phase ``bounds``.
        :param iso: iso-value (default ``0.0``).
        :param period: ``(Lx, Ly, Lz)`` if ``func`` is intrinsically periodic.
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
            field=func,
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
        :param scalars: name of the scalar array on the grid.
        :param iso: iso-value (default ``0.0``).
        :param name: phase name (auto-generated if omitted).
        """
        if scalars not in grid.point_data:
            err_msg = (
                f"Grid has no point scalar named {scalars!r}. "
                f"Available: {list(grid.point_data.keys())}"
            )
            raise ValueError(err_msg)

        nx, ny, nz = grid.dimensions
        pts = np.asarray(grid.points).reshape((nx, ny, nz, 3), order="F")
        xmin, ymin, zmin = pts[0, 0, 0]
        xmax, ymax, zmax = pts[-1, -1, -1]
        bounds = (
            float(xmin),
            float(xmax),
            float(ymin),
            float(ymax),
            float(zmin),
            float(zmax),
        )
        scalar = np.asarray(grid[scalars]).reshape((nx, ny, nz), order="F")
        dx = (xmax - xmin) / (nx - 1) if nx > 1 else 1.0
        dy = (ymax - ymin) / (ny - 1) if ny > 1 else 1.0
        dz = (zmax - zmin) / (nz - 1) if nz > 1 else 1.0

        def _nearest_field(
            x: np.ndarray,
            y: np.ndarray,
            z: np.ndarray,
        ) -> np.ndarray:
            xa = np.asarray(x)
            ya = np.asarray(y)
            za = np.asarray(z)
            ix = np.clip(np.round((xa - xmin) / dx).astype(int), 0, nx - 1)
            iy = np.clip(np.round((ya - ymin) / dy).astype(int), 0, ny - 1)
            iz = np.clip(np.round((za - zmin) / dz).astype(int), 0, nz - 1)
            return scalar[ix, iy, iz]

        instance = cls(
            field=_nearest_field,
            bounds=bounds,
            iso=iso,
            name=name,
            resolution=max(nx, ny, nz),
        )
        # Seed the grid cache so subsequent .grid() returns the original
        # (avoids resampling the field via nearest-neighbour).
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
            generate_binary_structure,
            label,
        )

        sg = self.grid()
        res = self._resolution
        scalar = np.asarray(sg[_IMPLICIT_SCALAR]).reshape((res, res, res), order="F")
        inside = scalar < self._iso
        # 26-connectivity (faces + edges + corners) — the physically correct
        # notion of "connected solid region" for a sampled continuous field.
        # scipy's default is 6-connectivity (faces only), which fragments
        # convex polyhedra at sharp vertices: a corner voxel and its diagonal
        # neighbors all sample slightly inside but only touch corner-to-corner,
        # so 6-connectivity reports them as separate components.  Genuinely
        # disjoint regions (e.g. periodic-split fragments) are separated by
        # many voxels and remain distinct under any connectivity.
        structure = generate_binary_structure(3, 3)
        labels, n_labels = label(inside, structure=structure)
        if n_labels == 0:
            return []

        xmin, xmax, ymin, ymax, zmin, zmax = self._bounds  # type: ignore[misc]
        dx = (xmax - xmin) / (res - 1)
        dy = (ymax - ymin) / (res - 1)
        dz = (zmax - zmin) / (res - 1)
        cell_volume = dx * dy * dz

        coms = center_of_mass(inside, labels=labels, index=range(1, n_labels + 1))
        slices = find_objects(labels)

        out: list[Piece] = []
        for label_id in range(1, n_labels + 1):
            mask = labels == label_id
            voxel_count = int(mask.sum())
            com_voxel = coms[label_id - 1]
            com_world = (
                xmin + com_voxel[0] * dx,
                ymin + com_voxel[1] * dy,
                zmin + com_voxel[2] * dz,
            )
            sl = slices[label_id - 1]
            piece_bounds = (
                xmin + sl[0].start * dx,
                xmin + (sl[0].stop - 1) * dx,
                ymin + sl[1].start * dy,
                ymin + (sl[1].stop - 1) * dy,
                zmin + sl[2].start * dz,
                zmin + (sl[2].stop - 1) * dz,
            )
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
            sg = self.grid()
            res = self._resolution
            scalar = np.asarray(sg[_IMPLICIT_SCALAR]).reshape(
                (res, res, res), order="F"
            )
            inside = scalar < self._iso
            if not inside.any():
                err_msg = "center_of_mass: field is positive everywhere on the grid"
                raise ValueError(err_msg)
            pts = np.asarray(sg.points).reshape((res, res, res, 3), order="F")
            mask = inside.astype(np.float64)
            denom = float(mask.sum())
            com = (
                float((mask * pts[..., 0]).sum() / denom),
                float((mask * pts[..., 1]).sum() / denom),
                float((mask * pts[..., 2]).sum() / denom),
            )
            return np.array(com)
        err_msg = "Cannot compute center_of_mass on an empty Phase"
        raise ValueError(err_msg)

    @cached_property
    def inertia_matrix(self: Phase) -> npt.NDArray[np.float64]:
        """Inertia tensor (about the origin) of the phase.

        Field-backed: grid quadrature.  CAD-backed: ``BRepGProp``.
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
            sg = self.grid()
            res = self._resolution
            scalar = np.asarray(sg[_IMPLICIT_SCALAR]).reshape(
                (res, res, res), order="F"
            )
            inside = scalar < self._iso
            if not inside.any():
                err_msg = "inertia_matrix: field is positive everywhere on the grid"
                raise ValueError(err_msg)
            pts = np.asarray(sg.points).reshape((res, res, res, 3), order="F")
            x = pts[..., 0][inside]
            y = pts[..., 1][inside]
            z = pts[..., 2][inside]
            xmin, xmax, ymin, ymax, zmin, zmax = self._bounds  # type: ignore[misc]
            dx = (xmax - xmin) / (res - 1)
            dy = (ymax - ymin) / (res - 1)
            dz = (zmax - zmin) / (res - 1)
            cell_volume = dx * dy * dz
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
        (:meth:`CadShape.copy` → ``BRepBuilderAPI_Copy``), the surface mesh
        and any cached grid are deep-copied, so an ``inplace=False``
        transform on the copy never touches the original.
        """
        new = Phase(
            field=self._field,
            bounds=self._bounds,
            iso=self._iso,
            period=self._period,
            name=self.name,
            resolution=self._resolution,
        )
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
    # ``inplace=False`` (default) returns a new Phase (operating on a
    # :meth:`copy`); ``inplace=True`` mutates ``self``.  Both return a
    # Phase.  The ``_apply_*`` workers mutate every live representation
    # (field, CAD, surface mesh, cached grid) via ``self`` so the result
    # stays coherent regardless of which mode the caller picked.
    # ------------------------------------------------------------------

    def translate(
        self: Phase, offset: Sequence[float], *, inplace: bool = False
    ) -> Phase:
        """Translate the phase by ``offset``.

        :param offset: ``(dx, dy, dz)`` world-space shift.
        :param inplace: mutate ``self`` when ``True``; otherwise (default)
            return a new translated :class:`Phase`.
        """
        if self.is_empty:
            err_msg = "Cannot translate an empty Phase"
            raise ValueError(err_msg)
        dx, dy, dz = float(offset[0]), float(offset[1]), float(offset[2])
        target = self if inplace else self.copy()
        target._apply_translate(dx, dy, dz)  # noqa: SLF001
        return target

    def _apply_translate(self: Phase, dx: float, dy: float, dz: float) -> None:
        if self._field is not None and self._bounds is not None:
            f = self._field

            def _translated_field(x, y, z, _f=f, _dx=dx, _dy=dy, _dz=dz):  # noqa: ANN001, ANN202
                return _f(x - _dx, y - _dy, z - _dz)

            b = self._bounds
            self._field = _translated_field
            self._bounds = (
                b[0] + dx,
                b[1] + dx,
                b[2] + dy,
                b[3] + dy,
                b[4] + dz,
                b[5] + dz,
            )
        if self._cad is not None:
            self._cad = self._cad.translate((dx, dy, dz))
        if self._surface_mesh is not None:
            self._surface_mesh = self._surface_mesh.translate(
                (dx, dy, dz), inplace=False
            )
        if self._cached_grid is not None:
            self._cached_grid = self._cached_grid.translate((dx, dy, dz), inplace=False)
        self._invalidate_derived()

    def scale(
        self: Phase,
        factor: float | tuple[float, float, float],
        *,
        inplace: bool = False,
    ) -> Phase:
        """Scale the phase about its center.

        Field-backed phases scale about their volumetric
        :attr:`center_of_mass`; CAD/mesh-backed phases about the BRep/mesh
        center.  Uniform ``factor`` is exact; a per-axis tuple rescales the
        bbox and (for an SDF) leaves the iso-distance approximate.

        :param factor: uniform scalar or ``(sx, sy, sz)``.
        :param inplace: mutate ``self`` when ``True``; otherwise (default)
            return a new scaled :class:`Phase`.
        """
        if self.is_empty:
            err_msg = "Cannot scale an empty Phase"
            raise ValueError(err_msg)
        if isinstance(factor, (int, float)):
            sx = sy = sz = float(factor)
        else:
            sx, sy, sz = (float(s) for s in factor)
        # One pivot shared by every representation so they stay coherent.
        pivot = self._scale_pivot()
        target = self if inplace else self.copy()
        target._apply_scale(sx, sy, sz, pivot)  # noqa: SLF001
        return target

    def _scale_pivot(self: Phase) -> tuple[float, float, float]:
        if self._field is not None:
            return tuple(float(v) for v in self.center_of_mass)  # type: ignore[return-value]
        if self._cad is not None:
            c = self._cad.center()
            return (float(c.x), float(c.y), float(c.z))
        c = self._surface_mesh.center  # type: ignore[union-attr]
        return (float(c[0]), float(c[1]), float(c[2]))

    def _apply_scale(
        self: Phase,
        sx: float,
        sy: float,
        sz: float,
        pivot: tuple[float, float, float],
    ) -> None:
        px, py, pz = pivot
        if self._field is not None and self._bounds is not None:
            f = self._field
            s_min = min(sx, sy, sz)

            def _scaled_field(  # noqa: ANN001, ANN202, PLR0913
                x, y, z, _f=f, _sx=sx, _sy=sy, _sz=sz, _px=px, _py=py, _pz=pz, _sm=s_min
            ):
                return (
                    _f(
                        (x - _px) / _sx + _px,
                        (y - _py) / _sy + _py,
                        (z - _pz) / _sz + _pz,
                    )
                    * _sm
                )

            b = self._bounds
            self._field = _scaled_field
            self._bounds = (
                px + (b[0] - px) * sx,
                px + (b[1] - px) * sx,
                py + (b[2] - py) * sy,
                py + (b[3] - py) * sy,
                pz + (b[4] - pz) * sz,
                pz + (b[5] - pz) * sz,
            )
        if self._cad is not None:
            from .cad import transform_geometry  # noqa: PLC0415

            matrix = np.array(
                [
                    [sx, 0.0, 0.0, px - sx * px],
                    [0.0, sy, 0.0, py - sy * py],
                    [0.0, 0.0, sz, pz - sz * pz],
                ],
                dtype=np.float64,
            )
            self._cad = transform_geometry(self._cad, matrix)
        if self._surface_mesh is not None:
            self._surface_mesh = self._surface_mesh.scale(
                (sx, sy, sz), point=(px, py, pz), inplace=False
            )
        # A cached input grid can't be rescaled in lockstep with the SDF
        # pivot; drop it so .grid() re-samples the scaled field.
        self._cached_grid = None
        self._invalidate_derived()

    def rotate(
        self: Phase,
        rotation: Rotation | np.ndarray,
        point: Sequence[float] | None = None,
        *,
        inplace: bool = False,
    ) -> Phase:
        """Rotate the phase by a SciPy :class:`~scipy.spatial.transform.Rotation`.

        - CAD-backed: delegates to ``CadShape.rotate`` (axis-angle form).
        - Field-backed: composes ``f'(p) = f(R^{-1}(p - point) + point)`` so
          the iso-surface rotates with the rest; ``bounds`` becomes the AABB
          of the rotated original AABB and ``period`` resets to ``None``.

        :param rotation: a ``Rotation`` object (or a 3x3 matrix), following
            PyVista's ``RotationLike``.  Build it with ``Rotation.from_euler``
            / ``from_matrix`` / ``from_quat`` / ``from_rotvec``.
        :param point: pivot point; defaults to the world origin.
        :param inplace: mutate ``self`` when ``True``; otherwise (default)
            return a new rotated :class:`Phase`.
        """
        if self.is_empty:
            err_msg = "Cannot rotate an empty Phase"
            raise ValueError(err_msg)
        rot = (
            rotation
            if isinstance(rotation, Rotation)
            else Rotation.from_matrix(np.asarray(rotation, dtype=np.float64))
        )
        pivot = (
            (0.0, 0.0, 0.0)
            if point is None
            else (float(point[0]), float(point[1]), float(point[2]))
        )
        target = self if inplace else self.copy()
        target._apply_rotate(rot, pivot)  # noqa: SLF001
        return target

    def _apply_rotate(
        self: Phase, rot: Rotation, pivot: tuple[float, float, float]
    ) -> None:
        px, py, pz = pivot
        p = np.array([px, py, pz], dtype=np.float64)
        if self._field is not None and self._bounds is not None:
            inv = rot.inv().as_matrix()
            rot_matrix = rot.as_matrix()
            f = self._field

            def _rotated_field(x, y, z, _f=f, _m=inv, _px=px, _py=py, _pz=pz):  # noqa: ANN001, ANN202
                xx = x - _px
                yy = y - _py
                zz = z - _pz
                lx = _m[0, 0] * xx + _m[0, 1] * yy + _m[0, 2] * zz + _px
                ly = _m[1, 0] * xx + _m[1, 1] * yy + _m[1, 2] * zz + _py
                lz = _m[2, 0] * xx + _m[2, 1] * yy + _m[2, 2] * zz + _pz
                return _f(lx, ly, lz)

            b = self._bounds
            corners = np.array(
                [
                    [x, y, z]
                    for x in (b[0], b[1])
                    for y in (b[2], b[3])
                    for z in (b[4], b[5])
                ]
            )
            rc = (corners - p) @ rot_matrix.T + p
            self._field = _rotated_field
            self._bounds = (
                float(rc[:, 0].min()),
                float(rc[:, 0].max()),
                float(rc[:, 1].min()),
                float(rc[:, 1].max()),
                float(rc[:, 2].min()),
                float(rc[:, 2].max()),
            )
            self._period = None  # rotation breaks axis-aligned periodicity
        if self._cad is not None:
            rotvec = rot.as_rotvec(degrees=True)
            angle = float(np.linalg.norm(rotvec))
            if angle != 0.0:
                axis = rotvec / angle
                self._cad = self._cad.rotate((px, py, pz), tuple(axis), angle)
        if self._surface_mesh is not None:
            self._surface_mesh = self._surface_mesh.rotate(
                rot, point=(px, py, pz), inplace=False
            )
        if self._cached_grid is not None:
            self._cached_grid = self._cached_grid.rotate(
                rot, point=(px, py, pz), inplace=False
            )
        self._invalidate_derived()

    def tile(self: Phase, rve: Rve, grid: tuple[int, int, int]) -> Phase:
        """Periodically tile the phase on the RVE over ``grid`` copies.

        A **producer**: always returns a new :class:`Phase` (no ``inplace``).
        Builds ``∏ grid`` translated copies of the current phase and
        fuses them.  Only implemented for CAD-backed phases today (the
        field-backed equivalent — domain folding via ``mod`` — lands in
        a follow-up; the periodic-shape work in
        :mod:`microgen.shape.implicit_ops` already covers it for
        :class:`Shape` directly).
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

        center = np.array(self._cad.center().to_tuple())
        copies = []
        for idx in np.ndindex(*grid):
            pos = center - rve.dim * (0.5 * np.array(grid) - 0.5 - np.array(idx))
            copies.append(translate_solid(self._cad.wrapped, pos))
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
