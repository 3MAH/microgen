from abc import ABC, abstractmethod
from typing import TYPE_CHECKING
import numpy.typing as npt
import numpy as np
from scipy.spatial.transform import Rotation
from microgen.shape import Shape
from microgen import (
    Rve,
    Phase,
    Cylinder,
    Box,
    fuse_shapes,
    periodic,
)

if TYPE_CHECKING:
    import cadquery as cq
    import pyvista as pv

    from microgen.shape import Vector3DType

class AbstractLattice(Shape):
    """
    Abstract Class to create strut-based lattice
    """

    def __init__(self,
                 strut_radius: float = 0.05,
                 cell_size: float = 1.0,
                 **kwargs: Vector3DType,
                 ) -> None:
        """
        Abstract Class to create strut-based lattice.
        The lattice will be created in a cube which size can be modified with 'cell_size'.
        The number of repetitions in each direction of the created geometry can be modified with 'repeat_cell'.

        :param center: center of the lattice
        :param orientation: orientation of the lattice
        :param strut_radius: radius of the struts
        :param cell_size: size of the cubic rve in which the lattice cell is enclosed
        
        """
        super.__init__(**kwargs)

        self.strut_radius = strut_radius
        self.cell_size = cell_size

        self.rve = Rve(dim=self.cell_size, center=self.center)

        self.vertices = self._compute_vertices()
        self.strut_centers = self._compute_strut_centers()
        self.strut_directions_cartesian = self._compute_strut_directions()
        self.strut_directions_euler = self._compute_euler_angles()

    @property
    @abstractmethod
    def strut_number(self) -> int: ...

    @property
    @abstractmethod
    def strut_height(self) -> float: ...

    @abstractmethod
    def _compute_vertices(self) -> npt.NDArray[np.float64]: ...

    @abstractmethod
    def _compute_strut_centers(self) -> npt.NDArray[np.float64]: ...

    @abstractmethod
    def _compute_strut_directions(self) -> npt.NDArray[np.float64]: ...

    def _compute_euler_angles(self) -> npt.NDArray[np.float64]:
        """Computes euler angles from default (1.0, 0.0, 0.0) oriented cylinder for all struts in the lattice"""

        default_direction = np.array([1.0, 0.0, 0.0])

        rotation_vector_array = np.zeros((self.strut_number, 3))
        euler_angles_array = np.zeros((self.strut_number, 3))

        for i in range(self.strut_number):
            if np.all(self.strut_directions_cartesian[i] == default_direction) or np.all(self.strut_directions_cartesian[i] == -default_direction):
                euler_angles_array[i] = np.zeros(3)
            else:
                axis = np.cross(default_direction, self.strut_directions_cartesian[i])
                axis /= np.linalg.norm(axis)
                angle = np.arccos(np.dot(default_direction, self.strut_directions_cartesian[i]))
                rotation_vector_array[i] = angle * axis
                euler_angles_array[i] = Rotation.from_rotvec(rotation_vector_array[i]).as_euler('zxz', degrees=True)

        return euler_angles_array

    def generate(self) -> cq.Shape:
        list_phases : list[Phase] = []
        list_periodic_phases : list[Phase] = []

        for i in range(self.strut_number):
            strut = Cylinder(
                center=tuple(self.strut_centers[i]),
                orientation=(self.strut_directions_euler[i, 2], self.strut_directions_euler[i, 1],
                             self.strut_directions_euler[i, 0]),
                height=self.strut_height,
                radius=self.strut_radius,
            )
            list_phases.append(Phase(strut.generate()))

        for phase_strut in list_phases:
            periodic_phase = periodic(phase=phase_strut, rve=self.rve)
            list_periodic_phases.append(periodic_phase)

        lattice = fuse_shapes([phase.shape for phase in list_periodic_phases], retain_edges=False)

        bounding_box = Box(center=self.center, dim_x=self.cell_size, dim_y=self.cell_size, dim_z=self.cell_size).generate()

        cut_lattice = bounding_box.intersect(lattice)

        return cut_lattice

    @property
    def volume(self) -> float:
        volume = self.generate().Volume()

        return volume