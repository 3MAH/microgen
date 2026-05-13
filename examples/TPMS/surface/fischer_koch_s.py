from pathlib import Path

from microgen import Tpms
from microgen.shape.surface_functions import fischer_koch_s

geometry = Tpms(surface_function=fischer_koch_s).with_repeat_cell(5).with_offset(0)
mesh = geometry.generate_surface_mesh(type_part="surface")

vtk_file = Path(__file__).parent / "surface.vtk"
mesh.save(vtk_file)

# pl = pv.Plotter()
# pl.add_mesh(mesh)
# pl.save_graphic('fischer_koch_s.pdf')
