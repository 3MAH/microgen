import numpy as np

from microgen import Phase, Rve, rasterPhase
from microgen.shape import Sphere


def test_phase_sphere_rasterize_must_have_correct_number_of_solids() -> None:
    """Rasterize Sphere by 3x3x3 grid must have 27 solids"""
    rve = Rve(dim=1)
    sphere = Sphere(radius=0.5).generate()
    grid = [3 for _ in range(3)]

    phase = Phase(shape=sphere)
    phases = phase.rasterize(rve=rve, grid=grid)
    raster = rasterPhase(phase=phase, rve=rve, grid=grid, phasePerRaster=False)

    assert len(phases) == np.prod(grid)
    assert len(phases) == len(raster.solids)
