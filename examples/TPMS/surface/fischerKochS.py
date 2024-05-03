from pathlib import Path

from microgen import Tpms
from microgen.shape.surface_functions import fischerKochS

geometry = Tpms(surface_function=fischerKochS, repeat_cell=5, offset=0)
mesh = geometry.generate_vtk(type_part="surface")

vtk_file = Path(__file__).parent / "surface.vtk"
mesh.save(vtk_file)

# pl = pv.Plotter()
# pl.add_mesh(mesh)
# pl.save_graphic('fischerKochS.pdf')
