"""Tests for mesh module."""

from pathlib import Path

import cadquery as cq
import numpy as np
import pyvista as pv

from microgen import Phase, Rve, Sphere, mesh, rasterPhase

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/

MESH_QUALITY_THRESHOLD = 0.5


def test_mesh_rastered_sphere_must_have_correct_number_of_cells() -> None:
    """Rasterize Sphere by 3x3x3 grid must have 27 cells."""
    Path("tests/data").mkdir(parents=True, exist_ok=True)

    grid = [3 for _ in range(3)]
    phases = rasterPhase(
        phase=Phase(shape=Sphere(radius=0.5).generate()),
        rve=Rve(),
        grid=grid,
    )
    compound = cq.Compound.makeCompound(
        [solid for phase in phases for solid in phase.solids],
    )
    cq.exporters.export(compound, fname="tests/data/compound.step")

    mesh(
        mesh_file="tests/data/compound.step",
        listPhases=phases,
        size=0.03,
        order=1,
        output_file="tests/data/compound.vtk",
    )

    vtk_mesh = pv.read("tests/data/compound.vtk")
    assert vtk_mesh["CellEntityIds"].shape[0] == vtk_mesh.n_cells > 0
    assert len(np.unique(vtk_mesh["CellEntityIds"])) == np.prod(grid)

    mesh_quality = np.mean(vtk_mesh.compute_cell_quality()["CellQuality"])
    assert mesh_quality > MESH_QUALITY_THRESHOLD
