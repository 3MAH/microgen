import os
import subprocess
import sys

import meshio

import microgen

USE_MMG = False
try:
    subprocess.run(["mmg3d_O3", "-h"])
    USE_MMG = True
except:
    print("mmg command did not work, check if it is installed or contact a developer")

if USE_MMG:
    meshIni = "Mesh.msh"

    if "data" not in os.listdir("."):
        os.mkdir("data")

    mesh = meshio.read(meshIni)
    mesh.write("data/meshIni.mesh")

    microgen.external.Mmg.mmg3d(input="data/meshIni.mesh", output="data/intermesh.mesh")
    microgen.external.Mmg.mmg3d(
        input="data/intermesh.mesh", output="data/finalmesh.mesh", ls=True, hsiz=0.03
    )

    meshFinal = meshio.read("data/finalmesh.mesh")
    meshFinal.write("finalmesh.msh", file_format="gmsh22")
    meshFinal.write("finalmesh.vtk")
