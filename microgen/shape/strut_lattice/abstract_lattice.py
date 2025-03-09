from __future__ import annotations
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, List
import numpy.typing as npt
import numpy as np
from scipy.spatial.transform import Rotation
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
import cadquery as cq
import pyvista as pv
from microgen.shape import Shape
from microgen import (
    Rve,
    Phase,
    Cylinder,
    Box,
    fuse_shapes,
    periodic,
    rotate,
)

if TYPE_CHECKING:
    from microgen.shape import Vector3DType, KwargsGenerateType

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
        super().__init__(**kwargs)

        self.strut_radius = strut_radius
        self.cell_size = cell_size

        self.rve = Rve(dim=self.cell_size, center=self.center)

        self.vertices = self._compute_vertices()
        self.strut_centers = self._compute_strut_centers()
        self.strut_directions_cartesian = self._compute_strut_directions()
        self.strut_rotations = self._compute_rotations()

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

    def _compute_rotations(self) -> List[Rotation]:
        """Computes euler angles from default (1.0, 0.0, 0.0) oriented cylinder for all struts in the lattice"""

        default_direction = np.array([1.0, 0.0, 0.0])

        rotations_list = []

        for i in range(self.strut_number):
            if np.all(self.strut_directions_cartesian[i] == default_direction) or np.all(self.strut_directions_cartesian[i] == -default_direction):
                rotation_vector = np.zeros(3)
                rotations_list.append(Rotation.from_rotvec(rotation_vector))
            else:
                axis = np.cross(default_direction, self.strut_directions_cartesian[i])
                axis /= np.linalg.norm(axis)
                angle = np.arccos(np.dot(default_direction, self.strut_directions_cartesian[i]))
                rotation_vector = angle * axis
                rotations_list.append(Rotation.from_rotvec(rotation_vector))#.as_euler('zxz', degrees=True)

        return rotations_list

    def generate(self, **_: KwargsGenerateType) -> cq.Shape:
        """Generate a strut-based lattice CAD shape using the given parameters."""
        list_phases : list[Phase] = []
        list_periodic_phases : list[Phase] = []

        for i in range(self.strut_number):
            strut = Cylinder(
                center=tuple(self.strut_centers[i]),
                orientation=self.strut_rotations[i],
                height=self.strut_height,
                radius=self.strut_radius,
            )
            list_phases.append(Phase(strut.generate()))

        for phase_strut in list_phases:
            periodic_phase = periodic(phase=phase_strut, rve=self.rve)
            list_periodic_phases.append(periodic_phase)

        lattice = fuse_shapes([phase.shape for phase in list_periodic_phases], retain_edges=False)

        bounding_box = Box(center=self.center, dim_x=self.cell_size, dim_y=self.cell_size, dim_z=self.cell_size).generate()

        #cut_lattice = bounding_box.intersect(lattice)
        cut_lattice = BRepAlgoAPI_Section(bounding_box.wrapped, lattice.wrapped).Shape()

        return cq.Shape(cut_lattice)
    
    @property
    def volume(self) -> float:
        volume = self.generate().Volume()

        return volume
    
    def generate_vtk(self, resolution: int = 100, **_: KwargsGenerateType) -> pv.PolyData:
        """Generate a strut-based lattice VTK shape using the given parameters."""
        lattice_structure = None
        
        for i in range(self.strut_number):
            strut = rotate(pv.Cylinder(
                center=tuple(self.strut_centers[i]),
                direction=(1.0, 0.0, 0.0),
                radius=self.strut_radius,
                height=self.strut_height,
                resolution=resolution,
                capping=True,  
            ).triangulate(), center=self.strut_centers[i], rotation=self.strut_rotations[i])
            
            if lattice_structure is None:
                lattice_structure = strut
            else:
                lattice_structure.boolean_union(strut)
        
        bounding_box = Box(center=self.center, dim_x=self.cell_size, dim_y=self.cell_size, dim_z=self.cell_size).generate_vtk()
        
        cut_lattice = bounding_box.boolean_intersection(lattice_structure)
        
        return cut_lattice

    def generateVtk( # noqa: N802
        self,
        resolution: int = 100,
        **kwargs: KwargsGenerateType,
        )-> pv.PolyData:
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk(
            resolution=resolution,
            **kwargs,
        )