"""Generate a Gyroid with callable thickness grading using microgen[meshers]."""

from microgen import Tpms
from microgen.shape.surface_functions import gyroid


def thickness_offset(x, y, z):
    """Vary the offset from 0.9 to 1.3 over the unit cell's z extent."""
    return 1.1 + 0.4 * z


geometry = Tpms(surface_function=gyroid, offset=thickness_offset)
mesh = geometry.generate_volume_mesh(
    backend="meshers",
    cells=24,
    order=1,
    geometry_tolerance=0.002,
    periodic=(True, True, False),  # Thickness varies along z, so z is not periodic.
    threads=1,
)
mesh.save("graded_gyroid_meshers.vtu")
print(f"Graded Gyroid: {mesh.n_cells} tetrahedra, {mesh.n_points} nodes")
