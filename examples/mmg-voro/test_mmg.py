import os

import microgen
import meshio

meshIni = "Mesh.msh"

if "data" not in os.listdir("."):
    os.mkdir("data")

mesh = meshio.read(meshIni)
mesh.write("data/meshIni.mesh")

microgen.mmg.mmg3d(input="data/meshIni.mesh", output="data/intermesh.mesh")
microgen.mmg.mmg3d(input="data/intermesh.mesh", output="data/finalmesh.mesh", ls=True, hgrad=1.1, hsiz=0.02) #, hgrad=1.1, hsiz=0.01)

meshFinal = meshio.read("data/finalmesh.mesh")
meshFinal.write("finalmesh.msh", file_format="gmsh22")
