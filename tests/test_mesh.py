import os

import cadquery as cq
import numpy as np
import pyvista as pv

from microgen import Phase, Rve, Sphere, mesh, rasterPhase


def test_mesh_rastered_sphere_must_have_correct_number_of_cells():
    """Rasterize Sphere by 3x3x3 grid must have 27 cells"""
    os.makedirs("tests/data", exist_ok=True)  # if data folder doesn't exist yet

    grid = [3 for _ in range(3)]
    phases = rasterPhase(
        phase=Phase(shape=Sphere(radius=0.5).generate()),
        rve=Rve(),
        grid=grid,
    )
    compound = cq.Compound.makeCompound(
        [solid for phase in phases for solid in phase.solids]
    )
    cq.exporters.export(compound, fname="tests/data/compound.step")

    mesh(
        mesh_file="tests/data/compound.step",
        list_phases=phases,
        size=0.03,
        order=1,
        output_file="tests/data/compound.vtk",
    )

    vtk_mesh = pv.read("tests/data/compound.vtk")
    assert vtk_mesh["CellEntityIds"].shape[0] == vtk_mesh.n_cells > 0
    assert len(np.unique(vtk_mesh["CellEntityIds"])) == np.prod(grid)
    assert np.mean(vtk_mesh.compute_cell_quality()["CellQuality"]) > 0.5
