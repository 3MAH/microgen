from pathlib import Path

import cadquery as cq
import numpy as np

from microgen import Box, ExtrudedPolygon, Phase, cutPhaseByShapeList, mesh

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

box = Box(dim_x=thickness, dim_y=60, dim_z=60)

shapeList = []
for seed in seedList:
    poly = ExtrudedPolygon(
        center=(seed[0] - thickness, seed[1], seed[2]),
        listCorners=[
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

honeycomb = cutPhaseByShapeList(phaseToCut=boxPhase, cqShapeList=shapeList)

step_file = str(Path(__file__).parent / "honeycomb.step")
stl_file = str(Path(__file__).parent / "honeycomb.stl")
cq.exporters.export(honeycomb.shape, step_file)
cq.exporters.export(honeycomb.shape, stl_file)
vtk_file = str(Path(__file__).parent / "honeycomb.vtk")
mesh(
    mesh_file=step_file,
    listPhases=[honeycomb],
    size=1,
    order=1,
    output_file=vtk_file,
)
