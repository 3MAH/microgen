from microgen import Tpms, tpms
import pyvista as pv

geometry = Tpms(
    surface_function=tpms.fischerKochS,
    repeat_cell=5
)
mesh = geometry.generateSurfaceVtk()

pl = pv.Plotter()
pl.add_mesh(mesh)
pl.save_graphic('fischerKochS.pdf')
