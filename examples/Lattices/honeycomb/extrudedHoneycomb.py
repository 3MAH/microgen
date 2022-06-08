from microgen import Box, ExtrudedPolygon, Phase, cutPhaseByShapeList
import numpy as np
import cadquery as cq

side_length = 2.5  # side in mm of the hexagon
poly_height = 2.5  # height in mm of the hexagon
theta = 30 * np.pi / 180  # half angle of the hexagone
density = 0.5  # relative density (roh*/roh_s) of the honeycomb

h0 = 0.5 * poly_height
h1 = np.cos(theta) * side_length
h2 = abs(np.sin(theta) * side_length)

# t = density*(2*h1*(h/l + abs(np.sin(theta))))/(h/l + 2)

thickness = 30  # mm


g = open("seedList.data", 'r', encoding="iso-8859-15")  # open data file

seedList = [[1, 1, 1]]


seedList = np.genfromtxt(g, delimiter='\t')

box = Box(dim_x=thickness, dim_y=60, dim_z=60)

shapeList = []
for seed in seedList:
    poly = ExtrudedPolygon(center=(seed[0] - thickness, seed[1], seed[2]),
                           listCorners=[(0, h2 + h0), (h1, h0),
                                        (h1, -h0), (0, -h2 - h0),
                                        (-h1, -h0), (-h1, h0),
                                        (0, h2 + h0)],
                           height=thickness)
    shapeList.append(poly.generate())

# generate CAD geometry
denseSample = Phase(shape=box.generate())

sample = cutPhaseByShapeList(phaseToCut=denseSample, cqShapeList=shapeList)

cq.exporters.export(sample.shape, 'poly1.step')
cq.exporters.export(sample.shape, 'honeycomb.stl')
# mesh(mesh_file='poly1.step', listPhases=[sample.getSolids()], size=0.1, order=1, output_mesh='Mesh.msh')
