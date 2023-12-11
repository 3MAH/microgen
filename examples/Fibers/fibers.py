import os

import cadquery as cq
import numpy as np
import progressbar  # pip install progressbar2

from microgen import Cylinder


def read_csv(filename):
    x = []
    y = []
    z = []
    diameter = []
    phi = []
    theta = []
    aspect_ratio = []
    index = []
    with open(filename) as f:
        f.readline()
        f.readline()
        for line in f:
            data = line.split(",")
            x.append(float(data[0]))
            y.append(float(data[1]))
            z.append(float(data[2]))
            diameter.append(float(data[3]))
            phi.append(float(data[4]))
            theta.append(float(data[5]))
            aspect_ratio.append(float(data[8]))
            index.append(int(data[9]))
    return x, y, z, diameter, phi, theta, aspect_ratio, index


def convert_angles(phi, theta):
    new_psi = phi
    new_theta = [90 for _ in phi]
    new_phi = theta

    return new_psi, new_theta, new_phi


x, y, z, diameter, phi, theta, aspect_ratio, index = read_csv("fibers.csv")
psi, theta, phi = convert_angles(phi, theta)


shapes = []
bar = progressbar.ProgressBar(max_value=len(x))
for i in range(len(x)):
    height = 250.0  # constant length = 150-350 µm (250) ou aléatoire
    radius = 0.5 * diameter[i]
    elem = Cylinder(
        center=(x[i], y[i], z[i]),
        orientation=(psi[i], theta[i], phi[i]),
        height=height,
        radius=radius,
    )
    shapes.append(elem.generate())
    bar.update(i)

assemb = cq.Assembly()
for shape in shapes:
    assemb.add(shape)


compound = assemb.toCompound()

cq.exporters.export(compound, "compound.step")
cq.exporters.export(compound, "compound.stl")


# mesh(mesh_file='compound.step', listPhases=raster[1], size=0.03, order=1, output_file='Mesh.msh')
