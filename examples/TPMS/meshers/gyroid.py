"""Generate a periodic Gyroid volume mesh with microgen[meshers]."""

from pathlib import Path

from microgen import Tpms
from microgen.shape.surface_functions import gyroid

geometry = Tpms(surface_function=gyroid, offset=1.1)
mesh = geometry.generate_volume_mesh(
    backend="meshers",
    cells=24,
    order=1,
    geometry_tolerance=0.002,
    periodic=(True, True, True),
    threads=1,
)
mesh.save(Path(__file__).parent / "gyroid_meshers.vtu")
print(f"Gyroid: {mesh.n_cells} tetrahedra, {mesh.n_points} nodes")
