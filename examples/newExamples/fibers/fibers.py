import os
import numpy as np
import cadquery as cq
from microgen import removeEmptyLines, Rve, BasicGeometry, periodic, cutParts, mesh, meshPeriodic, rasterShapeList, fuseParts


def read_csv(filename):
    x = []
    y = []
    z = []
    diameter = []
    phi = []
    theta = []
    aspect_ratio = []
    index = []
    with open(filename, 'r') as f:
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


listPhases = []
n = len(x) # 100
for i in range(0, n):
    # if i==4:
    #     continue
    # height = aspect_ratio[i]*diameter[i]
    height = 250. # constant length = 150-350 µm (250) ou aléatoire
    radius = 0.5*diameter[i]
    elem = BasicGeometry(number=index[i], shape="cylinder",
                        xc=x[i], yc=y[i], zc=z[i],
                        psi=psi[i], theta=theta[i], phi=phi[i],
                        param_geom={"height": height, "radius": radius}, 
                        path_data='data')
    listPhases.append(elem.generate())

assemb = cq.Assembly()
for phase in listPhases:
    assemb.add(phase)


compound = assemb.toCompound()

cq.exporters.export(compound, 'compound.step')
cq.exporters.export(compound, 'compound.stl')



# mesh(mesh_file='compound.step', listPhases=raster[1], size=size_mesh, order=1)