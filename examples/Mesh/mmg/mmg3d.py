import os
from pathlib import Path

import meshio

import microgen

data_dir = Path(__file__).parent / "data"
os.makedirs(data_dir, exist_ok=True)

msh_file = str(Path(__file__).parent / "Mesh.msh")

mesh = meshio.read(msh_file)
mesh.write(str(data_dir / "meshIni.mesh"))

microgen.external.Mmg.mmg3d(
    input=str(data_dir / "meshIni.mesh"),
    output=str(data_dir / "intermesh.mesh"),
)
microgen.external.Mmg.mmg3d(
    input=str(data_dir / "intermesh.mesh"),
    output=str(data_dir / "finalmesh.mesh"),
    ls=True,
    hsiz=0.03,
)

meshFinal = meshio.read(str(data_dir / "finalmesh.mesh"))
final_msh = str(Path(__file__).parent / "finalmesh.msh")
final_vtk = str(Path(__file__).parent / "finalmesh.vtk")
meshFinal.write(final_msh, file_format="gmsh22")
meshFinal.write(final_vtk)
