import numpy as np
import cadquery as cq

from microgen.Functions import removeEmptyLines
from microgen.Box import Box
from microgen.Sphere import Sphere
from microgen.Cylinder import Cylinder
from microgen.ExtrudedPolygon import ExtrudedPolygon
from microgen.Ellipsoid import Ellipsoid
from microgen.Capsule import Capsule
from microgen.Tpms import Tpms
from microgen.Polyhedron import Polyhedron


class BasicGeometry:
    def __init__(
        self, number, shape, xc, yc, zc, psi, theta, phi, param_geom, path_data=None
    ):
        self.number = number
        self.shape = shape
        self.xc = xc
        self.yc = yc
        self.zc = zc
        self.psi = psi
        self.theta = theta
        self.phi = phi
        self.param_geom = param_geom
        self.path_data = path_data

        self.center = np.array([self.xc, self.yc, self.zc])
        self.angle = np.array([self.psi, self.theta, self.phi])
        self.name = self.shape + str(self.number)

        if self.shape.lower() == "box":
            self.geometry = Box(
                center=self.center,
                angle=self.angle,
                a1=self.param_geom["a1"],
                a2=self.param_geom["a2"],
                a3=self.param_geom["a3"],
                number=self.number,
            )
        if self.shape.lower() == "cylinder":
            self.geometry = Cylinder(
                center=self.center,
                angle=self.angle,
                height=self.param_geom["height"],
                radius=self.param_geom["radius"],
                number=self.number,
            )
        if self.shape.lower() == "extrudedpolygon":
            self.geometry = ExtrudedPolygon(
                center=self.center,
                angle=self.angle,
                listCorners=self.param_geom["listCorners"],
                height=self.param_geom["height"],
                number=self.number,
            )
        if self.shape.lower() == "capsule":
            self.geometry = Capsule(
                center=self.center,
                angle=self.angle,
                height=self.param_geom["height"],
                radius=self.param_geom["radius"],
                number=self.number,
            )
        if self.shape.lower() == "sphere":
            self.geometry = Sphere(
                center=self.center,
                radius=self.param_geom["radius"],
                number=self.number
            )
        if self.shape.lower() == "ellipsoid":
            self.geometry = Ellipsoid(
                center=self.center,
                angle=self.angle,
                a1=self.param_geom["a1"],
                a2=self.param_geom["a2"],
                a3=self.param_geom["a3"],
                number=self.number,
            )
        if self.shape.lower() == "tpms":
            if self.param_geom["type_surface"] == "custom":
                self.geometry = Tpms(
                    center=self.center,
                    angle=self.angle,
                    type_surface=self.param_geom["type_surface"],
                    type_part=self.param_geom["type_part"],
                    thickness=self.param_geom["thickness"],
                    number=self.number,
                    function=self.param_geom["function"]
                )
            else:
                self.geometry = Tpms(
                    center=self.center,
                    angle=self.angle,
                    type_surface=self.param_geom["type_surface"],
                    type_part=self.param_geom["type_part"],
                    thickness=self.param_geom["thickness"],
                    number=self.number
                )

        if self.shape.lower() == "polyhedron":
            self.geometry = Polyhedron(
                dic=self.param_geom,
                number=self.number
            )

    def __cmp__(self, other):
        # return cmp(self.number, other.number)
        return (self.number > other.number) - (self.number < other.number)  # replacement for cmp function not availbale with Python3

    # ----------GENERATE PHASES----------------------------------------------------------------------------------

    def generate(self, rve=None):

        if self.shape.lower() == "box":
            cqshape = self.geometry.createBox()
        elif self.shape.lower() == "cylinder":
            cqshape = self.geometry.createCylinder()
        elif self.shape.lower() == "extrudedpolygon":
            cqshape = self.geometry.createExtrudedpolygon()
        elif self.shape.lower() == "capsule":
            cqshape = self.geometry.createCapsule()
        elif self.shape.lower() == "sphere":
            cqshape = self.geometry.createSphere()
        elif self.shape.lower() == "ellipsoid":
            cqshape = self.geometry.createEllipsoid()
        elif self.shape.lower() == "tpms":
            cqshape = self.geometry.createTpms(self.path_data, rve)
        elif self.shape.lower() == "polyhedron":
            cqshape = self.geometry.createPolyhedron()
        else:
            raise ValueError(self.shape + " is not recognised")

        return cq.Shape(cqshape.val().wrapped)


# def readPhases(path_data, phases_file, phases):

#     nphases = 0
#     cnt_phase = 0
#     nprops = []
#     # buf = ""
#     path_inputfile = path_data + phases_file
#     removeEmptyLines(path_inputfile)

#     try:
#         fp = open(path_inputfile)
#         for line in enumerate(fp):
#             if line == "\n":
#                 nphases -= 1
#             else:
#                 row_split = line[1].split()
#                 if nphases > 0:

#                     nprops_phase = 0
#                     if row_split[1] == "matrix":
#                         nprops_phase = 0
#                     elif row_split[1] == "sphere":
#                         nprops_phase = 1
#                     elif row_split[1] == "cylinder":
#                         nprops_phase = 2
#                     elif row_split[1] == "bar":
#                         nprops_phase = 2
#                     elif row_split[1] == "ellipsoid":
#                         nprops_phase = 3
#                     elif row_split[1] == "tpms":
#                         nprops_phase = 3

#                     nprops.append(nprops_phase)
#                 nphases += 1
#     finally:
#         fp.seek(0)
#         fp.close()

#     print(path_data)
#     print(phases_file)

#     try:
#         fp = open(path_inputfile)
#         for line in enumerate(fp):
#             print(line)
#             if line == "\n":
#                 cnt_phase -= 1
#             else:
#                 row_split = line[1].split()
#                 if cnt_phase > 0:

#                     number = int(row_split[0])
#                     shape = row_split[1]
#                     xc = float(row_split[2])
#                     yc = float(row_split[3])
#                     zc = float(row_split[4])
#                     psi = float(row_split[5])
#                     theta = float(row_split[6])
#                     phi = float(row_split[7])

#                     props = []
#                     if shape == "tpms":
#                         for prop in range(0, nprops[cnt_phase - 1]):
#                             props.append(row_split[8 + prop])
#                     else:
#                         for prop in range(0, nprops[cnt_phase - 1]):
#                             props.append(float(row_split[8 + prop]))

#                     pha = phase(
#                         number, shape, xc, yc, zc, psi, theta, phi, props, path_data
#                     )
#                     phases.append(pha)
#                     print(phases[-1].shape)
#                 cnt_phase += 1
#     finally:
#         fp.close()
