from microgen import *

l = 2.5 # side in mm of the hexagon
h = 2.5 # height in mm of the hexagon
theta = 30*np.pi/180 # half angle of the hexagone 
density = 0.5 # relative density (roh*/roh_s) of the honeycomb

sampleThickness = 30 # mm

h1 = np.cos(theta)*l
h2 = abs(np.sin(theta)*l)

# t = density*(2*h1*(h/l + abs(np.sin(theta))))/(h/l + 2)




sampleThickness = 30


g = open("seedList.data", 'r', encoding="iso-8859-15") # open data file

seedList = [[1,1,1]]


seedList = np.genfromtxt(g, delimiter = '\t')

box = BasicGeometry(number=0,
                    shape='box',
                    xc=0, yc=0, zc=0,
                    psi=0, theta=0, phi=0,
                    param_geom={"dim_x": sampleThickness,
                                "dim_y": 60,
                                "dim_z": 60})

phases = []
for seed in seedList:
    poly = BasicGeometry(number=0,
                         shape='extrudedpolygon',
                         xc=seed[0] - sampleThickness, yc=seed[1], zc=seed[2],
                         psi=0, theta=0, phi=0,
                         param_geom={"listCorners": [(0, h2 + h/2), (h1, h/2), 
                                                     (h1, -h/2), (0, -h2 - h/2), 
                                                     (-h1, -h/2), (-h1, h/2), 
                                                     (0, h2 + h/2)],
                                     "height": sampleThickness})
    phases.append(poly.generate())
    

# generate CAD geometry
denseSample = box.generate()

sample = cutPhaseByShapeList(denseSample, phases)

cq.exporters.export(sample[0], 'poly1.step')
# mesh(mesh_file='poly1.step', listPhases=[sample[1]], size=0.1, order=1)

