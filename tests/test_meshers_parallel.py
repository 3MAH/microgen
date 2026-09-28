"""CPU budgets and real concurrent generation through the microgen API."""

from threading import Barrier
from types import SimpleNamespace

import numpy as np
import pytest

from microgen import Tpms, _meshers, generate_meshers_parallel, parallel
from microgen.shape.surface_functions import gyroid


def test_native_generation_defaults_to_automatic_workers(monkeypatch):
    received = []

    def generate(*args, **options):
        received.append(options["threads"])
        return SimpleNamespace(diagnostics={})

    monkeypatch.setattr(_meshers.meshers, "generate", generate)
    monkeypatch.setattr(_meshers.meshers, "generate_intersection", generate)
    _meshers.generate(lambda x, y, z: x, (-1, 1) * 3, 8)
    _meshers.generate(lambda x, y, z: x, (-1, 1) * 3, 8, threads=2)
    _meshers.generate({"field": lambda x, y, z: x}, (-1, 1) * 3, 8)
    assert received == [None, 2, None]


def test_batch_runs_concurrently_and_caps_nested_workers(monkeypatch):
    monkeypatch.setattr(parallel, "_available_threads", lambda: 4)
    rendezvous = Barrier(2, timeout=5)

    class Shape:
        def __init__(self, number):
            self.number = number

        def generate_meshers(self, *, threads, **options):
            assert threads == 2
            rendezvous.wait()
            return SimpleNamespace(number=self.number, diagnostics={})

    results = generate_meshers_parallel([Shape(i) for i in range(4)], threads=2)
    assert [r.number for r in results] == list(range(4))
    assert all(r.diagnostics["batch_workers"] == 2 for r in results)


def test_batch_divides_automatic_native_budget(monkeypatch):
    monkeypatch.setattr(parallel, "_available_threads", lambda: 8)

    class Shape:
        def generate_meshers(self, *, threads, **options):
            assert threads == 4
            return SimpleNamespace(diagnostics={})

    results = generate_meshers_parallel([Shape(), Shape()], threads=None)
    assert all(r.diagnostics["batch_workers"] == 2 for r in results)


def test_batch_rejects_shared_mutable_objects_and_invalid_budgets(monkeypatch):
    monkeypatch.setattr(parallel, "_available_threads", lambda: 4)
    shape = object()
    with pytest.raises(ValueError, match="distinct"):
        generate_meshers_parallel([shape, shape])
    for kwargs in [
        {"max_workers": 0},
        {"threads": -1},
        {"threads": 5},
        {"threads": True},
    ]:
        with pytest.raises(ValueError):
            generate_meshers_parallel([shape], **kwargs)
    assert generate_meshers_parallel([]) == []


def test_batch_propagates_generation_failure():
    class Shape:
        def generate_meshers(self, **options):
            raise RuntimeError("generation failed")

    with pytest.raises(RuntimeError, match="generation failed"):
        generate_meshers_parallel([Shape()])


def test_real_batch_preserves_meshes_and_periodicity(monkeypatch):
    monkeypatch.setattr(parallel, "_available_threads", lambda: 2)
    shapes = [Tpms(gyroid, offset=0.6, resolution=16) for _ in range(2)]
    options = {"periodic": (True,) * 3}
    reference = shapes[0].generate_meshers(threads=1, **options)
    results = generate_meshers_parallel(shapes, **options)
    for result in results:
        np.testing.assert_array_equal(result.points, reference.points)
        np.testing.assert_array_equal(result.tetrahedra, reference.tetrahedra)
        for actual, expected in zip(
            result.periodic_pairs, reference.periodic_pairs, strict=True
        ):
            np.testing.assert_array_equal(actual, expected)
        assert result.diagnostics["minimum_mmg_quality"] >= 0.1
        assert result.diagnostics["batch_native_threads"] == 1
