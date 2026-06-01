from pathlib import Path

from microgen import Neper, Phase, mesh
from microgen.cad import make_compound

# Import the Polyhedra from a Neper tessellation file, build a CAD compound,
# export STEP, and tetrahedralise via gmsh.
tess_file = str(Path(__file__).parent / "test1.tess")
polyhedra = Neper.voronoi_from_tess_file(tess_file)

shapes = [poly.generate_cad() for poly in polyhedra]

compound = make_compound(shapes)
step_file = str(Path(__file__).parent / "compound.step")
compound.export_step(step_file)

phases = [Phase.from_cad(shape) for shape in shapes]

vtk_file = str(Path(__file__).parent / "Voronoi.vtk")
mesh(
    mesh_file=step_file,
    list_phases=phases,
    size=0.05,
    order=1,
    output_file=vtk_file,
)
