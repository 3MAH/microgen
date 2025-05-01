import os.path

import pyvista as pv


def show_example(
    filename,
    screenshot=False,
    color=None,
    show_edges=False,
    screen_name=None,
    cmap="coolwarm",
    rotate=False,
    zoom=None,
):
    if not os.path.exists(filename):
        print(filename + " not found")
        return

    mesh = pv.read(filename)
    plotter = pv.Plotter(off_screen=screenshot)
    if rotate:
        plotter.camera_position = "zy"
        plotter.camera.position = (-1, 0.5, 0.5)
        plotter.camera.focal_point = (mesh.center[0], mesh.center[1], mesh.center[2])

    if zoom:
        plotter.camera.zoom(zoom)

    plotter.add_mesh(mesh, show_edges=show_edges, color=color, cmap=cmap)
    try:
        plotter.remove_scalar_bar()
    except Exception:
        pass

    if screenshot:
        if screen_name is None:
            screen_name = filename.split("/")[-1].split(".")[0] + ".png"
        screen_name = "examples/" + screen_name
        if os.path.exists(screen_name):
            print(screen_name + " already exists")
        plotter.screenshot(screen_name, transparent_background=True)
    else:
        plotter.show()
        print(plotter.camera.position)


example_path = "../../examples/"
show_example(filename=example_path + "BasicShapes/shapes/shapes.stl", screenshot=True)
show_example(filename=example_path + "BasicShapes/platon/platon.stl", screenshot=True)
show_example(
    filename=example_path + "Lattices/honeycomb/honeycomb.stl", screenshot=True
)
show_example(
    filename=example_path + "Lattices/octetTruss/octettruss.stl", screenshot=True
)
show_example(filename=example_path + "TPMS/tpmsShell/tpms_shell.stl", screenshot=True)
show_example(
    filename=example_path + "TPMS/tpmsSphere/tpms_sphere.stl",
    screenshot=True,
    # zoom=0.2
)
show_example(filename=example_path + "TPMS/gyroid/gyroid.stl", screenshot=True)
show_example(
    filename=example_path + "3Doperations/repeatGeom/repeated_geometry.stl",
    screenshot=True,
)
show_example(filename=example_path + "Lattices/lattices.stl", screenshot=True)
show_example(
    filename=example_path + "Lattices/auxetic_custom_lattice.stl", screenshot=True
)
show_example(
    filename=example_path + "3Doperations/rasterEllipsoid/rasterEllipsoid.vtk",
    screenshot=True,
    show_edges=True,
    rotate=True,
    zoom=0.75,
)
show_example(
    filename=example_path + "3Doperations/voronoi/Voronoi.vtk",
    screenshot=True,
    show_edges=False,
)
show_example(
    filename=example_path + "3Doperations/voronoiGyroid/Gyroid-voro.vtk",
    screenshot=True,
    show_edges=False,
)
show_example(
    filename=example_path + "3Doperations/voronoi/Voronoi.vtk",
    screenshot=True,
    show_edges=True,
    screen_name="Mesh.png",
)
show_example(
    filename=example_path + "Lattices/octetTruss/octettruss.vtk",
    screenshot=True,
    show_edges=True,
    color="white",
    cmap=None,
    screen_name="meshPeriodic.png",
)
show_example(
    filename=example_path + "Mesh/mmg/finalmesh.vtk",
    screenshot=True,
    show_edges=True,
    color="white",
    cmap=None,
    screen_name="mmg.png",
)
show_example(
    filename=example_path + "Mesh/remesh/data/initial_gyroid_mesh.vtk",
    screenshot=True,
    show_edges=True,
    color="white",
    cmap=None,
    screen_name="initial_mesh.png",
)
show_example(
    filename=example_path + "Mesh/remesh/data/remeshed_gyroid_mesh.vtk",
    screenshot=True,
    show_edges=True,
    color="white",
    cmap=None,
    screen_name="remeshed_mesh.png",
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
