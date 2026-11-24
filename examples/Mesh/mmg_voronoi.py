"""
Level-set remeshing of a Voronoi mesh
=====================================

Remesh a Voronoi mesh with ``mmg3d``, then discretize a level set
stored in the mesh, with a gradation of 1.1.

Requires the ``mmg`` executables: ``conda install -c conda-forge mmg``.
"""

from pathlib import Path

import meshio

import microgen

data_dir = Path("mmg_voronoi_data")
data_dir.mkdir(exist_ok=True)

msh_file = "mmg_voronoi.msh"

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
    hgrad=1.1,
    hsiz=0.02,
)

mesh_final = meshio.read(str(data_dir / "finalmesh.mesh"))
mesh_final.write(str(data_dir / "finalmesh.msh"), file_format="gmsh22")
mesh_final.write(str(data_dir / "finalmesh.vtk"))

# %%
# .. image:: /_static/examples/mmg-voro.png
#
# sphinx_gallery_thumbnail_path = "_static/examples/mmg-voro.png"
