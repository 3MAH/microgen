"""Tests for the Spinodoid shape."""

from __future__ import annotations

import numpy as np
import pyvista as pv
import pytest

from microgen import Spinodoid

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/
# ruff: noqa: E501 line-too-long


# ---------------------------------------------------------------------------
# Constructor sanity / basic shape
# ---------------------------------------------------------------------------


def test_spinodoid_constructor_default_fft_grf() -> None:
    """Default constructor builds a periodic GRF Spinodoid with target density."""
    sp = Spinodoid(k0=15.0, bandwidth=3.0, resolution=16, density=0.5, seed=0)
    assert sp._frep is not None
    assert sp._frep.modes.shape[1] == 3
    assert sp._frep.modes.dtype == np.int64
    assert len(sp._frep.modes) > 10
    assert isinstance(sp.grid, pv.StructuredGrid)
    # Porosity within Monte Carlo + clip noise of target
    cell_volume = float(np.prod(sp.cell_size * sp.repeat_cell))
    assert abs(abs(sp.grid_solid.volume) / cell_volume - 0.5) < 0.05


def test_spinodoid_invalid_density_raises() -> None:
    with pytest.raises(ValueError, match="density"):
        Spinodoid(density=0.0, resolution=12, seed=0)
    with pytest.raises(ValueError, match="density"):
        Spinodoid(density=1.5, resolution=12, seed=0)


def test_spinodoid_offset_density_mutually_exclusive() -> None:
    with pytest.raises(ValueError, match="cannot be given at the same time"):
        Spinodoid(offset=0.1, density=0.5, resolution=12, seed=0)
    with pytest.raises(ValueError, match="must be given"):
        Spinodoid(resolution=12, seed=0)


def test_spinodoid_offset_used_directly() -> None:
    """Passing `offset` skips the Monte Carlo search and uses it as the iso-value."""
    sp = Spinodoid(k0=15.0, bandwidth=3.0, resolution=16, offset=0.05, seed=0)
    assert sp._frep.threshold == 0.05


# ---------------------------------------------------------------------------
# F-rep callable + periodicity invariants
# ---------------------------------------------------------------------------


def test_spinodoid_func_evaluable_and_sign_matches_solid() -> None:
    """`_func` is evaluable; sign convention matches grid_solid membership."""
    sp = Spinodoid(k0=15.0, bandwidth=3.0, resolution=16, density=0.5, seed=0)
    rng = np.random.default_rng(0)
    pts = rng.uniform(0.0, sp.cell_size[0], size=(100, 3))
    values = sp._func(pts[:, 0], pts[:, 1], pts[:, 2])
    assert values.shape == (100,)
    # Negative ↔ inside the solid (Shape F-rep convention)
    inside_count = int((values < 0).sum())
    # Allow ±15% slack for Monte Carlo noise on a 100-sample test
    assert 35 <= inside_count <= 65


def test_spinodoid_field_is_bit_exact_periodic() -> None:
    """`_func(x + L) == _func(x)` to numerical noise."""
    sp = Spinodoid(k0=15.0, bandwidth=3.0, resolution=16, density=0.5, seed=0)
    Lx, Ly, Lz = sp.cell_size.tolist()
    rng = np.random.default_rng(1)
    x = rng.uniform(0.0, Lx, size=10)
    y = rng.uniform(0.0, Ly, size=10)
    z = rng.uniform(0.0, Lz, size=10)
    v0 = sp._func(x, y, z)
    vx = sp._func(x + Lx, y, z)
    vy = sp._func(x, y + Ly, z)
    vz = sp._func(x, y, z + Lz)
    assert np.max(np.abs(vx - v0)) < 1e-10
    assert np.max(np.abs(vy - v0)) < 1e-10
    assert np.max(np.abs(vz - v0)) < 1e-10


def test_spinodoid_seed_reproducibility() -> None:
    """Same seed → identical mode arrays; different seeds → different."""
    sp_a = Spinodoid(k0=15.0, bandwidth=3.0, resolution=12, density=0.5, seed=42)
    sp_b = Spinodoid(k0=15.0, bandwidth=3.0, resolution=12, density=0.5, seed=42)
    sp_c = Spinodoid(k0=15.0, bandwidth=3.0, resolution=12, density=0.5, seed=99)
    np.testing.assert_array_equal(sp_a._frep.modes, sp_b._frep.modes)
    np.testing.assert_array_equal(sp_a._frep.amplitudes, sp_b._frep.amplitudes)
    assert not np.array_equal(sp_a._frep.modes, sp_c._frep.modes) or len(
        sp_a._frep.modes
    ) != len(sp_c._frep.modes)


# ---------------------------------------------------------------------------
# generate() — CAD shape
# ---------------------------------------------------------------------------


def test_spinodoid_generate_returns_solid_or_compound() -> None:
    """generate() builds a Solid or Compound with volume matching grid_solid."""
    pytest.importorskip("OCP")
    from OCP.TopAbs import TopAbs_COMPOUND, TopAbs_SOLID

    sp = Spinodoid(k0=15.0, bandwidth=3.0, resolution=16, density=0.5, seed=0)
    shape = sp.generate_cad()
    shape_type = shape.wrapped.ShapeType()
    assert shape_type in (TopAbs_SOLID, TopAbs_COMPOUND)
    # Volume should match grid_solid within 5% (mesh→CAD roundtrip noise)
    grid_vol = abs(sp.grid_solid.volume)
    cad_vol = abs(shape.volume())
    assert abs(cad_vol - grid_vol) / max(grid_vol, 1e-12) < 0.05


def test_spinodoid_generate_step_export_is_3d() -> None:
    """STEP export imports back as a 3D entity in gmsh."""
    pytest.importorskip("OCP")
    pytest.importorskip("gmsh")
    import tempfile

    import gmsh

    sp = Spinodoid(k0=15.0, bandwidth=3.0, resolution=12, density=0.5, seed=0)
    shape = sp.generate_cad()
    with tempfile.NamedTemporaryFile(suffix=".step", delete=False) as f:
        step_path = f.name
    shape.export_step(step_path)

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 0)
    try:
        gmsh.model.occ.importShapes(step_path, highestDimOnly=True)
        gmsh.model.occ.synchronize()
        n3d = len(gmsh.model.getEntities(3))
    finally:
        gmsh.finalize()
    assert n3d >= 1


# ---------------------------------------------------------------------------
# Sampling grid
# ---------------------------------------------------------------------------


def test_spinodoid_grid_has_cap_endpoints() -> None:
    """Structured grid extrema land bit-exactly at cell boundaries."""
    sp = Spinodoid(k0=15.0, bandwidth=3.0, resolution=16, density=0.5, seed=0)
    bounds = sp.grid.bounds
    Lx, Ly, Lz = sp.cell_size.tolist()
    assert bounds[0] == 0.0
    assert bounds[1] == Lx
    assert bounds[2] == 0.0
    assert bounds[3] == Ly
    assert bounds[4] == 0.0
    assert bounds[5] == Lz


def test_spinodoid_repeat_cell_tiling() -> None:
    """repeat_cell extends the grid bbox to repeat_cell × cell_size."""
    sp = Spinodoid(
        k0=15.0,
        bandwidth=3.0,
        resolution=8,
        density=0.5,
        cell_size=1.0,
        repeat_cell=2,
        seed=0,
    )
    bounds = sp.grid.bounds
    assert bounds[1] == 2.0
    assert bounds[3] == 2.0
    assert bounds[5] == 2.0


# ---------------------------------------------------------------------------
# F-rep integration with other Shapes
# ---------------------------------------------------------------------------


def test_spinodoid_frep_boolean_op_with_implicit_sphere() -> None:
    """Spinodoid is a first-class F-rep Shape: boolean ops produce a new Shape."""
    from microgen.shape import Shape, from_field

    sp = Spinodoid(k0=15.0, bandwidth=3.0, resolution=12, density=0.5, seed=0)

    def sphere_sdf(
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
    ) -> np.ndarray:
        return np.sqrt((x - 0.5) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2) - 0.3

    sphere = from_field(sphere_sdf, bounds=(0.0, 1.0, 0.0, 1.0, 0.0, 1.0))
    intersected = sp & sphere
    assert isinstance(intersected, Shape)
    # Smoke test: callable evaluable on arbitrary points.
    pts = np.array([[0.5, 0.5, 0.5]])
    val = intersected._func(pts[:, 0], pts[:, 1], pts[:, 2])
    assert val.shape == (1,)
