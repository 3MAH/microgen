import microgen
import meshio

meshIni = "Mesh.msh"

mesh = meshio.read(meshIni)
mesh.write("meshIni.mesh")

microgen.mmg.mmg3d(input="meshIni.mesh", output="intermesh.mesh")
microgen.mmg.mmg3d(input="intermesh.mesh", output="finalmesh.mesh", ls=True, hsiz=0.03)

meshFinal = meshio.read("finalmesh.mesh")
meshFinal.write("finalmesh.msh", file_format="gmsh22")
