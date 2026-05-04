from pathlib import Path

import numpy as np

from microgen import Box, ExtrudedPolygon, Phase, cut_phase_by_shape_list, mesh
from microgen.mesh import MeshOptions

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
    shapeList.append(poly.generate())

boxPhase = Phase(shape=box.generate())

honeycomb = cut_phase_by_shape_list(phase_to_cut=boxPhase, shapes=shapeList)

step_file = str(Path(__file__).parent / "honeycomb.step")
stl_file = str(Path(__file__).parent / "honeycomb.stl")
honeycomb.shape.export_step(step_file)
honeycomb.shape.export_stl(stl_file)
vtk_file = str(Path(__file__).parent / "honeycomb.vtk")
mesh(
    mesh_file=step_file,
    list_phases=[honeycomb],
    options=MeshOptions(size=1, order=1, output_file=vtk_file),
)
