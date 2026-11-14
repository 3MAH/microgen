"""Sampled offsets, a cylindrical sector and infill with microgen[meshers]."""

import numpy as np
import pyvista as pv

from microgen import CylindricalTpms, Infill, Tpms
from microgen.shape.surface_functions import gyroid

sampled = Tpms(gyroid, offset=1.1, resolution=16)
sampled.offset = 1.1 + 0.2 * sampled.grid.points[:, 2]
sampled.generate_volume_mesh(backend="meshers", cells=16).save("sampled_gyroid.vtu")

# One angular cell: opposite sector faces are related by a rotation.
sector = CylindricalTpms(radius=2, surface_function=gyroid, offset=1.1, resolution=8)
sector_mesh = sector.generate_volume_mesh(
    backend="meshers", cells=16, periodic=(False, True, False)
)
sector_mesh.save("gyroid_sector.vtu")
print("Sector rotation/translation matrix:")
print(sector_mesh.field_data["periodic_transform_y"])

envelope = pv.Sphere(radius=0.45, theta_resolution=32, phi_resolution=32)
infill = Infill(
    obj=envelope, surface_function=gyroid, offset=1.1, repeat_cell=1, resolution=12
)
infill.generate_volume_mesh(backend="meshers", cells=24).save("gyroid_infill.vtu")
