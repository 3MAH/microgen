"""The MMG selector is shorthand for the existing VTK/remesh sequence."""

import importlib
import shutil
import sys

import numpy as np
import pytest

from microgen import Sphere, Spinodoid, Tpms
from microgen.shape.surface_functions import gyroid


@pytest.mark.parametrize("kind", ["shape", "tpms", "spinodoid"])
def test_mmg_backend_forwards_the_existing_vtk_mesh_and_options(monkeypatch, kind):
    shapes = {
        "shape": lambda: Sphere(radius=0.35, center=(2, 3, 4)),
        "tpms": lambda: Tpms(gyroid, offset=1.1, resolution=8),
        "spinodoid": lambda: Spinodoid(offset=0, resolution=8, seed=4),
    }
    shape = shapes[kind]()
    seed_options = {"resolution": 8} if kind == "shape" else {}
    if kind == "tpms":
        seed_options = {"type_part": "lower skeletal", "algo_resolution": 8}
    expected = shape.generate_volume_mesh(**seed_options)
    result = expected.copy()
    calls = []

    def remesh(seed, **options):
        np.testing.assert_array_equal(seed.points, expected.points)
        np.testing.assert_array_equal(seed.cells, expected.cells)
        calls.append(options)
        return result

    monkeypatch.setitem(sys.modules, "meshers", None)
    monkeypatch.setattr(
        importlib.import_module("microgen.remesh"),
        "remesh_keeping_boundaries_for_fem",
        remesh,
    )
    options = {"periodic": False, "hausd": 0.002, "hmax": 0.1, "tol": 1e-7}
    actual = shape.generate_volume_mesh(backend="mmg", **seed_options, **options)
    assert actual is result and calls == [options]


def test_mmg_backend_retains_remesher_defaults_and_failures(monkeypatch):
    calls = []

    def remesh(seed, **options):
        calls.append(options)
        raise RuntimeError("existing MMG failure")

    monkeypatch.setattr(
        importlib.import_module("microgen.remesh"),
        "remesh_keeping_boundaries_for_fem",
        remesh,
    )
    with pytest.raises(RuntimeError, match="existing MMG failure"):
        Sphere(radius=0.35).generate_volume_mesh(backend="mmg", resolution=8)
    assert calls == [{}]


@pytest.mark.skipif(shutil.which("mmg3d_O3") is None, reason="MMG CLI not installed")
def test_mmg_backend_runs_the_real_pipeline():
    shape = Tpms(gyroid, offset=1.1, resolution=8)
    result = shape.generate_volume_mesh(backend="mmg", periodic=False, hmax=0.2)
    assert result.n_cells > 0 and result.volume > 0
    assert np.isfinite(result.points).all()
    assert set(result.cells_dict) == {10}  # VTK_TETRA
