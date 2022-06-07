import pyvista as pv
import os.path

stl_files = [
    "../../examples/newExamples/shapes/shapes.stl",
    "../../examples/newExamples/tpms_shell/tpms_shell.stl",
    "../../examples/newExamples/tpms_sphere/tpms_sphere.stl",
    "../../examples/octetTruss/octettruss.stl",
    "../../examples/honeycomb/honeycomb.stl",
    "../../examples/newExamples/tpms/tpms.stl",
    "../../examples/newExamples/repeatGeom/repeated_geometry.stl",
]

for filename in stl_files:
    image_name = filename.split('/')[-1].split('.')[0] + '.png'
    if not os.path.exists(filename):
        print(filename + " not found")
        continue
    elif os.path.exists(image_name):
        print(image_name + " already exists")
        continue
    else:
        mesh = pv.read(filename)
        plotter = pv.Plotter(off_screen=True)
        plotter.add_mesh(mesh)
        plotter.screenshot(image_name, transparent_background=True)