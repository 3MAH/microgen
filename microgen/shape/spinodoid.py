"""Spinodoid.

=============================================
Spinodoid (:mod:`microgen.shape.spinodoid`)
=============================================

A Spinodoid is a periodic stochastic microstructure built as a Gaussian
Random Field generated via FFT of white noise filtered by a Gaussian
shell. The Fourier coefficients of the filtered spectrum are kept as a
sparse expansion on the reciprocal lattice, so the implicit field is
continuously evaluable everywhere and participates in implicit boolean
ops on the same footing as TPMS shapes.

The implicit field is exactly periodic on ``cell_size × repeat_cell`` —
every kept mode lives on the reciprocal lattice (integer index in each
axis), so periodicity is a data-structure invariant, not a runtime flag.
"""

from __future__ import annotations

import contextlib
from collections.abc import Sequence
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import pyvista as pv

from ..operations import rotate
from ._frep_grf import _FrepGRF, _normalize_cell_size, compute_threshold_for_porosity
from .shape import Shape

if TYPE_CHECKING:
    from microgen.cad import CadShape
    from microgen.shape import KwargsGenerateType, Vector3DType


class Spinodoid(Shape):
    """Class to generate a periodic Spinodoid microstructure.

    Construction is in two stages:

    1. **Spectral**: an FFT of filtered white noise yields a sparse set of
       reciprocal-lattice modes ``(modes, amplitudes, phases)`` whose
       magnitudes survive a relative-amplitude threshold.
    2. **Geometric**: the iso-surface ``f(x) = offset`` defines the solid
       (``f(x) ≥ offset`` is inside). Either ``offset`` or ``density`` must
       be given, mirroring :class:`microgen.Tpms`:

       - ``offset``: the iso-value directly (a scalar in field units).
       - ``density``: target solid volume fraction in (0, 1); microgen
         finds the ``offset`` matching this fraction via Monte Carlo on
         the F-rep callable.

       The structured grid and ``grid_solid`` / ``surface`` properties are
       derived by sampling the callable on a ``(resolution+1)^3`` grid
       covering the full unit cell.

    The implicit field ``f(x) = Σᵢ Aᵢ cos(kᵢ·x + φᵢ) − offset`` is
    *bit-for-bit periodic* on ``cell_size`` (every ``kᵢ`` lives on the
    reciprocal lattice). Periodicity is an invariant — there is no
    ``periodic=True`` flag.

    :param k0: dominant wave number (Gaussian shell center).
    :param bandwidth: Gaussian shell width.
    :param offset: iso-value defining the solid surface (mutually
        exclusive with ``density``).
    :param density: target solid volume fraction in (0, 1) (mutually
        exclusive with ``offset``).
    :param anisotropy: ``(aₓ, a_y, a_z)`` axis-stretch factors applied
        to wavenumbers in the shell filter ``K_eff² = (aₓ Kₓ)² + …``.
    :param cell_size: unit-cell side length(s); scalar or 3-tuple.
    :param repeat_cell: integer tiling along each axis.
    :param resolution: FFT grid resolution per cell (the structured grid
        will have ``resolution + 1`` points per axis to cover both
        endpoints of the periodic cell exactly).
    :param mode_amplitude_threshold: keep FFT modes with
        ``|F̂| > threshold * max(|F̂|)``.
    :param seed: RNG seed; ``None`` for non-deterministic behavior.
    :param center: geometry center after construction (translation applied
        in :meth:`generate_surface_mesh` / :meth:`generate_cad`).
    :param orientation: rotation angles (Euler) applied at the end.
    """

    def __init__(  # noqa: PLR0913
        self: Spinodoid,
        k0: float = 30.0,
        bandwidth: float = 5.0,
        offset: float | None = None,
        density: float | None = None,
        anisotropy: Sequence[float] = (1.0, 1.0, 1.0),
        cell_size: float | Sequence[float] = 1.0,
        repeat_cell: int | Sequence[int] = 1,
        resolution: int = 32,
        mode_amplitude_threshold: float = 1e-3,
        seed: int | None = None,
        center: Vector3DType = (0.0, 0.0, 0.0),
        orientation: Vector3DType = (0.0, 0.0, 0.0),
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the Spinodoid object — see class docstring."""
        super().__init__(center=center, orientation=orientation, **kwargs)

        if offset is not None and density is not None:
            err_msg = (
                "offset and density cannot be given at the same time. Give only one."
            )
            raise ValueError(err_msg)
        if offset is None and density is None:
            err_msg = "offset or density must be given. Give one of them."
            raise ValueError(err_msg)
        if density is not None and not 0.0 < density < 1.0:
            err_msg = f"density must be in (0, 1). Given: {density}"
            raise ValueError(err_msg)

        self.k0 = float(k0)
        self.bandwidth = float(bandwidth)
        self.anisotropy = np.asarray(anisotropy, dtype=float)
        self.offset = float(offset) if offset is not None else None
        self.density = float(density) if density is not None else None
        self.resolution = int(resolution)
        self.mode_amplitude_threshold = float(mode_amplitude_threshold)
        self.seed = seed

        self._init_cell_parameters(cell_size, repeat_cell)

        self._frep: _FrepGRF
        self._setup_frep_field()

    def _init_cell_parameters(
        self: Spinodoid,
        cell_size: float | Sequence[float],
        repeat_cell: int | Sequence[int],
    ) -> None:
        self.cell_size = _normalize_cell_size(cell_size)

        if isinstance(repeat_cell, int):
            self.repeat_cell = np.array(
                [repeat_cell, repeat_cell, repeat_cell],
                dtype=int,
            )
        elif len(repeat_cell) == 3:
            self.repeat_cell = np.asarray(repeat_cell, dtype=int)
        else:
            err_msg = f"`repeat_cell` must have length 3. Given: {repeat_cell}"
            raise ValueError(err_msg)

    def _setup_frep_field(self: Spinodoid) -> None:
        """Build the F-rep modes, resolve the iso-value, wire ``_func``/``_bounds``."""
        self._frep = _FrepGRF.from_fft_grf(
            k0=self.k0,
            bandwidth=self.bandwidth,
            anisotropy=self.anisotropy,
            cell_size=self.cell_size,
            resolution=self.resolution,
            seed=self.seed,
            mode_amplitude_threshold=self.mode_amplitude_threshold,
        )

        if self.offset is not None:
            self._frep.threshold = self.offset
        else:
            self._frep.threshold = compute_threshold_for_porosity(
                self._frep,
                porosity=self.density,
                seed=self.seed,
            )

        # Shape uses "negative = inside"; _frep.evaluate is positive inside, so flip the sign.
        frep = self._frep

        def _signed_field(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            return -frep.evaluate(x, y, z)

        self._func = _signed_field
        lx, ly, lz = (self.cell_size * self.repeat_cell).tolist()
        self._bounds = (0.0, float(lx), 0.0, float(ly), 0.0, float(lz))

    @cached_property
    def grid(self: Spinodoid) -> pv.StructuredGrid:
        """Return the structured-grid sampling of the implicit field."""
        points, values = self._frep.sample_on_grid(
            resolution=self.resolution,
            repeat_cell=self.repeat_cell,
        )
        sg = pv.StructuredGrid(points[..., 0], points[..., 1], points[..., 2])
        sg["surface"] = values.ravel(order="F")
        return sg

    @cached_property
    def grid_solid(self: Spinodoid) -> pv.UnstructuredGrid:
        """Return the volumetric solid grid (cells where ``f >= 0``)."""
        return self.grid.clip_scalar(scalars="surface", value=0.0, invert=False)

    @cached_property
    def surface(self: Spinodoid) -> pv.PolyData:
        """Return the iso-surface mesh of the spinodoid."""
        return self.grid_solid.extract_surface().clean().triangulate()

    # ---- public generators -------------------------------------------------

    def generate_volume_mesh(
        self: Spinodoid,
        **_: KwargsGenerateType,
    ) -> pv.UnstructuredGrid:
        """Generate the volumetric solid grid (3D cells) with center+rotation applied."""
        grid_vol = self.grid_solid.copy()
        grid_vol = rotate(grid_vol, center=(0, 0, 0), rotation=self.orientation)
        return grid_vol.translate(xyz=self.center)

    def generate_surface_mesh(
        self: Spinodoid,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate the iso-surface as a PolyData with center+rotation applied."""
        polydata = self.surface.copy()
        polydata = rotate(polydata, center=(0, 0, 0), rotation=self.orientation)
        return polydata.translate(xyz=self.center)

    def generate_cad(
        self: Spinodoid,
        **_: KwargsGenerateType,
    ) -> CadShape:
        """Generate an OCCT CAD shape (Solid or Compound) of the spinodoid.

        Mirrors :meth:`microgen.Tpms.generate_cad`: extracts the structured-grid
        surface, builds a closed multi-face periodic shell with planar BREP
        faces on the cell sides plus per-triangle planar faces on the
        spinodoid interior, and applies the ``center`` / ``orientation``
        rigid transform. Spinodoids commonly decompose into multiple disjoint
        solid components — the return type may be a Compound of fused Solids.

        Requires the optional ``[cad]`` install extra.
        """
        mesh = self.surface
        pts = np.asarray(mesh.points, dtype=np.float64)
        tris = mesh.faces.reshape(-1, 4)[:, 1:].astype(np.int64)

        shape = _try_make_solid(_mesh_to_periodic_shell(pts, tris, self._bounds))
        shape = rotate(obj=shape, center=(0, 0, 0), rotation=self.orientation)
        shape = shape.translate(self.center)
        # Fallback for OCCT Volume() on solids it flags invalid (rigid transforms preserve volume).
        with contextlib.suppress(AttributeError, ValueError):
            shape._mesh_volume = float(abs(self.grid_solid.volume))
        return shape


def _mesh_to_periodic_shell(
    points: npt.NDArray[np.float64],
    triangles: npt.NDArray[np.int64],
    bounds: Sequence[float],
) -> CadShape:
    """Build a sewn OCCT shell from a triangulated surface inside an axis-aligned box.

    Triangles whose three vertices share the mesh's exact extremum on a cube
    plane are grouped per face and converted to a single planar BREP face via
    :func:`microgen.cad.mesh_to_planar_face`. Remaining triangles become one
    planar BREP face per triangle (raw, not pre-sewn — sewing must stitch
    them to cap wires along the seam). Everything is sewn into a closed shell.

    Mirrors :meth:`microgen.Tpms._mesh_to_periodic_shell` (tpms.py:683) but
    parameterised by ``bounds`` rather than an instance's
    ``cell_size × repeat_cell`` (Spinodoid lives at ``[0, L]^3``, Tpms at
    ``[-L/2, +L/2]^3``).
    """
    from OCP.BRepBuilderAPI import (
        BRepBuilderAPI_MakeEdge,
        BRepBuilderAPI_MakeFace,
        BRepBuilderAPI_MakeWire,
        BRepBuilderAPI_Sewing,
    )
    from OCP.gp import gp_Pnt
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer

    from microgen.cad import CadShape as _CadShape
    from microgen.cad import mesh_to_planar_face

    pts = np.asarray(points, dtype=np.float64)
    tris = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)

    drift_tol = 1e-9 * float(max(abs(b) for b in bounds) or 1.0)

    consumed = np.zeros(tris.shape[0], dtype=bool)
    on_plane: list[tuple[int, int, float, npt.NDArray[np.int64]]] = []
    for axis in range(3):
        for sign in (-1, +1):
            extremum = (
                float(pts[:, axis].max()) if sign > 0 else float(pts[:, axis].min())
            )
            expected = bounds[2 * axis + (1 if sign > 0 else 0)]
            if abs(extremum - expected) > drift_tol:
                continue
            vert_on = pts[:, axis] == extremum
            tri_on = np.all(vert_on[tris], axis=1) & ~consumed
            if tri_on.any():
                on_plane.append((axis, sign, extremum, np.where(tri_on)[0]))
                consumed |= tri_on
    interior_idx = np.where(~consumed)[0]

    bbox_diag = float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))
    sew_tol = max(1e-9, 1e-6 * bbox_diag)
    sewing = BRepBuilderAPI_Sewing(sew_tol)

    for axis, sign, extremum, tri_idx in on_plane:
        origin = [0.0, 0.0, 0.0]
        origin[axis] = extremum
        normal = [0.0, 0.0, 0.0]
        normal[axis] = float(sign)
        planar = mesh_to_planar_face(pts, tris[tri_idx], origin, normal)
        exp = TopExp_Explorer(planar.wrapped, TopAbs_FACE)
        while exp.More():
            sewing.Add(exp.Current())
            exp.Next()

    interior_tris = tris[interior_idx]
    used_vertices = np.unique(interior_tris)
    pnt_cache = {
        int(i): gp_Pnt(float(pts[i, 0]), float(pts[i, 1]), float(pts[i, 2]))
        for i in used_vertices
    }
    for a, b, c in interior_tris:
        ga, gb, gc = pnt_cache[int(a)], pnt_cache[int(b)], pnt_cache[int(c)]
        e1 = BRepBuilderAPI_MakeEdge(ga, gb).Edge()
        e2 = BRepBuilderAPI_MakeEdge(gb, gc).Edge()
        e3 = BRepBuilderAPI_MakeEdge(gc, ga).Edge()
        wire = BRepBuilderAPI_MakeWire(e1, e2, e3).Wire()
        face = BRepBuilderAPI_MakeFace(wire).Face()
        sewing.Add(face)

    sewing.Perform()
    return _CadShape(sewing.SewedShape())


def _try_make_solid(shape: CadShape) -> CadShape:
    """Best-effort upgrade of a sewn shell (or compound of shells) to a Solid.

    Mirrors :meth:`microgen.Tpms._try_make_solid` (tpms.py:1006). For each
    closed shell, builds a ``TopoDS_Solid`` via ``BRepBuilderAPI_MakeSolid``
    and reorients via ``BRepLib.OrientClosedSolid_s`` so the volume encloses
    the right region. Multiple disjoint shells fuse into a Compound.
    """
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeSolid
    from OCP.BRepLib import BRepLib
    from OCP.TopAbs import TopAbs_SHELL
    from OCP.TopExp import TopExp_Explorer

    from microgen.cad import CadShape as _CadShape
    from microgen.cad import _topods_cast
    from microgen.operations import fuse_shapes

    cast_shell = _topods_cast("Shell")
    cast_solid = _topods_cast("Solid")

    def _make_solid(shell) -> CadShape | None:
        try:
            solid = BRepBuilderAPI_MakeSolid(shell).Solid()
        except (ValueError, RuntimeError):
            return None
        BRepLib.OrientClosedSolid_s(cast_solid(solid))
        return _CadShape(solid)

    wrapped = shape.wrapped
    if wrapped.ShapeType() == TopAbs_SHELL:
        built = _make_solid(cast_shell(wrapped))
        return built if built is not None else shape

    exp = TopExp_Explorer(wrapped, TopAbs_SHELL)
    solids: list[CadShape] = []
    while exp.More():
        built = _make_solid(cast_shell(exp.Current()))
        if built is not None:
            solids.append(built)
        exp.Next()
    if not solids:
        return shape
    if len(solids) == 1:
        return solids[0]
    return fuse_shapes(solids, retain_edges=False)
