from microgen import CylindricalTpms
from microgen.shape.surface_functions import gyroid


def swapped_gyroid(x, y, z):
    return gyroid(x=z, y=y, z=x)


geometry = CylindricalTpms(
    surface_function=swapped_gyroid,
    offset=0.5,
    cell_size=(1, 1, 1),
    repeat_cell=(1, 0, 1),
    radius=1,
    resolution=20,
)
sheet = geometry.sheet

# plotter = pv.Plotter()
# plotter.add_mesh(sheet, color="w")
# plotter.view_xy()
# plotter.enable_parallel_projection()
# plotter.show_axes()
# plotter.show()
