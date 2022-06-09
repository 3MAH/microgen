import pyvista as pv
import os.path


def show_example(filename, screenshot=False, color=None, show_edges=False, screen_name=None, cmap='coolwarm'):
    if not os.path.exists(filename):
        print(filename + " not found")
        return

    mesh = pv.read(filename)
    plotter = pv.Plotter(off_screen=screenshot)
    plotter.add_mesh(mesh, show_edges=show_edges, color=color, cmap=cmap)
    try:
        plotter.remove_scalar_bar()
    except:
        pass

    if screenshot:
        if screen_name is None:
            screen_name = filename.split('/')[-1].split('.')[0] + '.png'
        screen_name = 'examples/' + screen_name
        if os.path.exists(screen_name):
            print(screen_name + " already exists")
        plotter.screenshot(screen_name, transparent_background=True)
    else:
        plotter.show()


example_path = "../../examples/"
show_example(
    filename=example_path + "BasicShapes/shapes/shapes.stl",
    screenshot=True
)
show_example(
    filename=example_path + "BasicShapes/platon/platon.stl",
    screenshot=True
)
show_example(
    filename=example_path + "Lattices/honeycomb/honeycomb.stl",
    screenshot=True
)
show_example(
    filename=example_path + "Lattices/octetTruss/octettruss.stl",
    screenshot=True
)
show_example(
    filename=example_path + "TPMS/tpms_shell/tpms_shell.stl",
    screenshot=True
)
show_example(
    filename=example_path + "TPMS/tpms_sphere/tpms_sphere.stl",
    screenshot=True
)
show_example(
    filename=example_path + "TPMS/gyroid/gyroid.stl",
    screenshot=True
)
show_example(
    filename=example_path + "3Doperations/repeatGeom/repeated_geometry.stl",
    screenshot=True
)
show_example(
    filename=example_path + "3Doperations/Voronoi/Voronoi.vtk",
    screenshot=True,
    show_edges=False
)
show_example(
    filename=example_path + "3Doperations/VoronoiGyroid/Gyroid-voro.vtk",
    screenshot=True,
    show_edges=False
)
show_example(
    filename=example_path + "3Doperations/Voronoi/Voronoi.vtk",
    screenshot=True,
    show_edges=True,
    screen_name="Mesh.png"
)
show_example(
    filename=example_path + "Lattices/octetTruss/octettruss.vtk",
    screenshot=True,
    show_edges=True,
    color='white',
    cmap=None,
    screen_name="meshPeriodic.png"
)


# show_example(
#     filename=example_path + "TPMS/tpms/tpms.stl",
#     screenshot=False
# )


stl_files = [
    "BasicShapes/shapes/shapes.stl",
    "BasicShapes/platon/platon.stl",
    "Lattices/honeycomb/honeycomb.stl",
    "Lattices/octettruss/octettruss.stl",
    "TPMS/tpms_shell/tpms_shell.stl",
    "TPMS/tpms_sphere/tpms_sphere.stl",
    "TPMS/tpms/tpms.stl",
    "TPMS/gyroid/gyroid.stl",
    "3Doperations/repeatGeom/repeated_geometry.stl",
    "3Doperations/Voronoi/Voronoi.vtk",
    "3Doperations/VoronoiGyroid/Gyroid-voro.vtk",
    # "3Doperations/rasterEllipsoid/rasterEllipsoid.vtk",
]
