from .abstract_lattice import AbstractLattice
from itertools import product
import numpy as np
import numpy.typing as npt
import math as m

_UNIT_CUBE_SIZE = 1.0
_NUM_CUBE_VERTICES = 8
_STRUT_NUMBER = 24
_STRUT_HEIGHTS = m.sqrt(2.0)/2.0

def _generate_base_vertices() -> npt.NDArray[np.float64]:
    cube_vertices = list(product([-_UNIT_CUBE_SIZE/2, _UNIT_CUBE_SIZE/2], repeat=3))
    
    face_centers = [
        [sign * 0.5 if i == axis else 0.0 for i in range(3)]
        for axis in range(3) for sign in [-1, 1]
    ]

    return np.array(cube_vertices + face_centers)

def _generate_face_to_vertices(vertices: np.ndarray) -> dict:
    """
    Dynamically generates a dictionary associating the indices of the face centers with the indices of the cube vertices.
    """
    cube_vertices = vertices[:_NUM_CUBE_VERTICES]
    face_centers = vertices[_NUM_CUBE_VERTICES:]

    face_to_vertices = {}

    for face_index, face_center in enumerate(face_centers):
        axis = np.argmax(np.abs(face_center))
        value = face_center[axis]

        connected_vertices = [
            i for i, vertex in enumerate(cube_vertices) if np.isclose(vertex[axis], value)
        ]

        face_to_vertices[face_index + _NUM_CUBE_VERTICES] = connected_vertices

    return face_to_vertices

def _generate_vertex_pairs(vertices: np.ndarray) -> np.ndarray:
    """
    Generates the pairs of indices of the fcc strut vertices between the face centers and the cube vertices.
    """
    face_to_vertices = _generate_face_to_vertices(vertices)

    vertex_pairs_indices = [
        [face_index, vertex_index]
        for face_index, vertex_indices in face_to_vertices.items()
        for vertex_index in vertex_indices
    ]

    return np.array(vertex_pairs_indices)

_VERTICES = _generate_base_vertices()
_STRUT_VERTEX_PAIRS = _generate_vertex_pairs(_VERTICES)

class FaceCenteredCubic(AbstractLattice):
    """
Class to create a unit face-centered cubic lattice of given cell size and density or strut radius
    """
    
    def __init__(self,
                 *args, **kwargs
                 ) -> None:

        super().__init__(*args, **kwargs, strut_number=_STRUT_NUMBER, strut_heights=_STRUT_HEIGHTS)

    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        """Compute the vertices of the face-centered cubic lattice."""
        return self.center + self.cell_size * _VERTICES

    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        """Compute the centers of the struts."""
        return np.mean(self.vertices[_STRUT_VERTEX_PAIRS], axis=1)

    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        """Compute the normalized direction vectors of the struts."""
        vectors = np.diff(self.vertices[_STRUT_VERTEX_PAIRS], axis=1).squeeze()
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)