"""
Level-set remeshing with mmg
============================

Remesh a gmsh mesh with ``mmg3d``, then discretize a level set
stored in the mesh.

Requires the ``mmg`` executables: ``conda install -c conda-forge mmg``.
"""

from pathlib import Path

import meshio

import microgen

data_dir = Path("mmg3d_data")
data_dir.mkdir(exist_ok=True)

msh_file = "mmg3d.msh"

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

mesh_final = meshio.read(str(data_dir / "finalmesh.mesh"))
mesh_final.write(str(data_dir / "finalmesh.msh"), file_format="gmsh22")
mesh_final.write(str(data_dir / "finalmesh.vtk"))

# %%
# .. image:: /_static/examples/mmg.png
#
# sphinx_gallery_thumbnail_path = "_static/examples/mmg.png"
