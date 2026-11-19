"""Tests for :func:`microgen.shape.periodic_shell.try_make_solid`."""

import numpy as np
import pytest
import pyvista as pv

# ruff: noqa: S101

pytest.importorskip("OCP")

from microgen.shape.periodic_shell import (  # noqa: E402
    mesh_to_periodic_shell,
    try_make_solid,
)

# Cell larger than the meshes: no triangle lies on a cap plane, so the
# shell is sewn from the mesh triangles alone.
_CELL = (-1.0, 4.0, -1.0, 2.0, -1.0, 2.0)


def _sewn(mesh: pv.PolyData):
    mesh = mesh.triangulate().clean()
    tris = mesh.faces.reshape(-1, 4)[:, 1:].astype(np.int64)
    return mesh_to_periodic_shell(
        np.asarray(mesh.points, dtype=np.float64), tris, _CELL
    )


def test_single_shell_becomes_solid_with_positive_volume() -> None:
    solid = try_make_solid(_sewn(pv.Box(bounds=(0, 1, 0, 1, 0, 1), level=1)))
    assert solid.volume() == pytest.approx(1.0, rel=1e-9)


def test_inward_normals_are_reoriented() -> None:
    box = pv.Box(bounds=(0, 1, 0, 1, 0, 1), level=1).triangulate()
    box = box.flip_faces() if hasattr(box, "flip_faces") else box.flip_normals()
    solid = try_make_solid(_sewn(box))
    assert solid.volume() == pytest.approx(1.0, rel=1e-9)


def test_disjoint_shells_become_compound_of_solids() -> None:
    two = pv.Box(bounds=(0, 1, 0, 1, 0, 1), level=1) + pv.Box(
        bounds=(2, 3, 0, 1, 0, 1), level=1
    )
    out = try_make_solid(_sewn(two))
    assert len(out.solids()) == 2
    assert out.volume() == pytest.approx(2.0, rel=1e-9)
