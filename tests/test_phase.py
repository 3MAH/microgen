from microgen.shape import Sphere
from microgen import Phase, Rve, rasterPhase

import numpy as np


def test_phase_sphere_rasterize():
    rve = Rve(dim_x=1, dim_y=1, dim_z=1)
    sphere = Sphere(radius=0.5).generate()
    grid = [3 for _ in range(3)]
    
    phases = Phase(shape=sphere).rasterize(rve=rve, grid=grid)

    assert len(phases) == np.prod(np.array(grid) - 1)