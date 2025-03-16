from .abstract_lattice import AbstractLattice
from itertools import product
import numpy as np
import numpy.typing as npt
import math as m

class FaceCenteredCubic(AbstractLattice):
    """
Class to create a unit face-centered cubic lattice of given cell size and density or strut radius
    """
    
    def __init__(self,
                 *args, **kwargs
                 ) -> None:

        self._base_vertices = self._generate_base_vertices()
        self._strut_vertex_pairs = self._generate_vertex_pairs()
        super().__init__(*args, **kwargs, strut_number=24, strut_heights=m.sqrt(2.0)/2.0)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        unit_cube_size = 1.0
        cube_vertices = list(product([-unit_cube_size/2, unit_cube_size/2], repeat=3))
        
        face_centers = [
            [sign * unit_cube_size/2 if i == axis else 0.0 for i in range(3)]
            for axis in range(3) for sign in [-1, 1]
        ]

        return np.array(cube_vertices + face_centers)

    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        """Compute the vertices of the face-centered cubic lattice."""
        base_vertices = self._generate_base_vertices()
        return self.center + self.cell_size * base_vertices
    
    def _generate_face_to_vertices(self) -> dict:
        """
        Dynamically generates a dictionary associating the indices of the face centers with the indices of the cube vertices.
        """
        num_cube_vertices = 8
        cube_vertices = self._base_vertices[:num_cube_vertices]
        face_centers = self._base_vertices[num_cube_vertices:]

        face_to_vertices = {}

        for face_index, face_center in enumerate(face_centers):
            axis = np.argmax(np.abs(face_center))
            value = face_center[axis]

            connected_vertices = [
                i for i, vertex in enumerate(cube_vertices) if np.isclose(vertex[axis], value)
            ]

            face_to_vertices[face_index + num_cube_vertices] = connected_vertices

        return face_to_vertices

    def _generate_vertex_pairs(self) -> np.ndarray:
        """
        Generates the pairs of indices of the fcc strut vertices between the face centers and the cube vertices.
        """
        face_to_vertices = self._generate_face_to_vertices()

        vertex_pairs_indices = [
            [face_index, vertex_index]
            for face_index, vertex_indices in face_to_vertices.items()
            for vertex_index in vertex_indices
        ]

        return np.array(vertex_pairs_indices)

    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        """Compute the centers of the struts."""
        return np.mean(self.vertices[self._strut_vertex_pairs], axis=1)

    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        """Compute the normalized direction vectors of the struts."""
        vectors = np.diff(self.vertices[self._strut_vertex_pairs], axis=1).squeeze()
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)