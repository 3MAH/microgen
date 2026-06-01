from pathlib import Path

import numpy as np

from microgen import Box, ExtrudedPolygon, Phase, cut_phase_by_shape_list, mesh

side_length = 2.5  # side in mm of the hexagon
poly_height = 2.5  # height in mm of the hexagon
theta = 30 * np.pi / 180  # half angle of the hexagone

h0 = 0.5 * poly_height
h1 = np.cos(theta) * side_length
h2 = abs(np.sin(theta) * side_length)

thickness = 30  # mm

data_file = str(Path(__file__).parent / "seedList.data")
with open(data_file) as f:
    seedList = [[1, 1, 1]]
    seedList = np.genfromtxt(f, delimiter="\t")

box = Box(dim=(thickness, 60, 60))

shapeList = []
for seed in seedList:
    poly = ExtrudedPolygon(
        center=(seed[0] - thickness, seed[1], seed[2]),
        list_corners=[
            (0, h2 + h0),
            (h1, h0),
            (h1, -h0),
            (0, -h2 - h0),
            (-h1, -h0),
            (-h1, h0),
            (0, h2 + h0),
        ],
        height=thickness,
    )
    shapeList.append(poly.generate_cad())

boxPhase = Phase.from_cad(box.generate_cad())

honeycomb = cut_phase_by_shape_list(phase_to_cut=boxPhase, shapes=shapeList)

step_file = str(Path(__file__).parent / "honeycomb.step")
stl_file = str(Path(__file__).parent / "honeycomb.stl")
honeycomb.cad.export_step(step_file)
honeycomb.cad.export_stl(stl_file)
vtk_file = str(Path(__file__).parent / "honeycomb.vtk")
mesh(
    mesh_file=step_file,
    list_phases=[honeycomb],
    size=1,
    order=1,
    output_file=vtk_file,
)
