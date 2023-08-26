from microgen import Tpms
from microgen.shape.surface_functions import fischerKochS
import pyvista as pv

geometry = Tpms(
    surface_function=fischerKochS,
    repeat_cell=5
)
mesh = geometry.generateVtk(type_part="surface")

mesh.save("surface.vtk")

# pl = pv.Plotter()
# pl.add_mesh(mesh)
# pl.save_graphic('fischerKochS.pdf')
