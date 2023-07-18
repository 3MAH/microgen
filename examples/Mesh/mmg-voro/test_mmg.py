import subprocess
import os

import microgen
import meshio

USE_MMG = False
try:
    subprocess.check_output("mmg3d_O3", stderr=subprocess.STDOUT)
    USE_MMG= True    
except(subprocess.CalledProcessError, FileNotFoundError):
    print(
        "mmg command did not work, check if it is installed or contact a developer"
    )
    USE_MMG= False 

if USE_MMG:
    meshIni = "Mesh.msh"

    if "data" not in os.listdir("."):
        os.mkdir("data")

    mesh = meshio.read(meshIni)
    mesh.write("data/meshIni.mesh")

    microgen.external.Mmg.mmg3d(input="data/meshIni.mesh", output="data/intermesh.mesh")
    microgen.external.Mmg.mmg3d(
        input="data/intermesh.mesh",
        output="data/finalmesh.mesh",
        ls=True,
        hgrad=1.1,
        hsiz=0.02,
    )

    meshFinal = meshio.read("data/finalmesh.mesh")
    meshFinal.write("finalmesh.msh", file_format="gmsh22")
    meshFinal.write("finalmesh.vtk")
