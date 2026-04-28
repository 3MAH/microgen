"""Strut-based lattice structures.

========================================
Strut Lattice (:mod:`microgen.shape.strut_lattice`)
========================================

.. jupyter-execute::
   :hide-code:

   import pyvista
   pyvista.set_jupyter_backend('static')
   pyvista.global_theme.background = 'white'
   pyvista.global_theme.window_size = [600, 400]
   pyvista.global_theme.axes.show = False
   pyvista.global_theme.smooth_shading = True
   pyvista.global_theme.split_sharp_edges = True

"""

from .abstract_lattice import AbstractLattice
from .body_centered_cubic import BodyCenteredCubic
from .cubic import Cubic
from .cuboctahedron import Cuboctahedron
from .custom_lattice import CustomLattice
from .diamond import Diamond
from .face_centered_cubic import FaceCenteredCubic
from .octahedron import Octahedron
from .octet_truss import OctetTruss
from .rhombic_cuboctahedron import RhombicCuboctahedron
from .rhombic_dodecahedron import RhombicDodecahedron
from .truncated_cube import TruncatedCube
from .truncated_cuboctahedron import TruncatedCuboctahedron
from .truncated_octahedron import TruncatedOctahedron

__all__ = [
    "AbstractLattice",
    "BodyCenteredCubic",
    "Cubic",
    "Cuboctahedron",
    "CustomLattice",
    "Diamond",
    "FaceCenteredCubic",
    "TruncatedCuboctahedron",
    "Octahedron",
    "OctetTruss",
    "TruncatedCube",
    "RhombicCuboctahedron",
    "RhombicDodecahedron",
    "TruncatedOctahedron",
]
