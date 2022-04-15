import cadquery as cq
from microgen import Rve, BasicGeometry

import argparse

surface   = "custom"
part      = "skeletal"
thickness = 0.1

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--func', nargs=1, 
                    help='TPMS custom function (e.g. python custom.py -f "cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)")')

args = parser.parse_args()

# function  = "cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)"
function = args.func[0]
print(function)


# Size of the mesh
size_mesh = 0.03
a = 1.0
b = 1.0
c = 1.0

revel = Rve(a, b, c, size_mesh)

elem = BasicGeometry(number=0, shape='tpms',
                     xc=0.5  , yc=0.5  , zc=0.5,
                     psi=0.  , theta=0., phi=0.,
                     param_geom={"type_surface": surface,
                                 "type_part"   : part,
                                 "thickness"   : thickness,
                                 "function"    : function},
                     path_data='data')


elem.geometry.createSurfaces(rve=revel,
                             sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03,
                             path_data='data')
part = elem.generate(rve=revel)

cq.exporters.export(part, 'part.step')