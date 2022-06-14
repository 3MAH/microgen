from microgen import Box, ExtrudedPolygon, Phase, cutPhaseByShapeList, mesh
import numpy as np
import cadquery as cq

side_length = 2.5  # side in mm of the hexagon
poly_height = 2.5  # height in mm of the hexagon
theta = 30 * np.pi / 180  # half angle of the hexagone

h0 = 0.5 * poly_height
h1 = np.cos(theta) * side_length
h2 = abs(np.sin(theta) * side_length)

thickness = 30  # mm

with open('seedList.data', 'r') as f:
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

cq.exporters.export(honeycomb.shape, "honeycomb.step")
cq.exporters.export(honeycomb.shape, "honeycomb.stl")
mesh(
    mesh_file="honeycomb.step",
    listPhases=[honeycomb],
    size=1,
    order=1,
    output_file="honeycomb.vtk",
)
