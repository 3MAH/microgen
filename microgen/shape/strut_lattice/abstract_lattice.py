from __future__ import annotations

from abc import abstractmethod
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import TYPE_CHECKING

import cadquery as cq
import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.spatial.transform import Rotation

from ...mesh import mesh, mesh_periodic
from ...operations import fuse_shapes
from ...periodic import periodic_split_and_translate
from ...phase import Phase
from ...rve import Rve
from ..box import Box
from ..cylinder import Cylinder
from ..shape import Shape
from ..sphere import Sphere

if TYPE_CHECKING:
    from microgen.shape import KwargsGenerateType, Vector3DType

# TODO add option to initialize lattice by giving density
# instead of strut radius


class AbstractLattice(Shape):
    """
    Abstract Class to create strut-based lattice
    """

    _UNIT_CUBE_SIZE = 1.0

    def __init__(
        self,
        strut_radius: float = 0.05,
        strut_heights: float | list[float] | None = None,
        base_vertices: npt.NDArray[np.float64] | None = None,
        strut_vertex_pairs: npt.NDArray[np.int64] | None = None,
        cell_size: float = 1.0,
        strut_joints: bool = False,
        **kwargs: Vector3DType | Rotation,
    ) -> None:
        """
        Abstract Class to create strut-based lattice.
        The lattice will be created in a cube which size can be
        modified with 'cell_size'.

        :param strut_radius: radius of the struts
        :param strut_height: either the unique height of all struts (float),
        or a list of strut heights (list[float]). Enter value for a size 1 rve.
        :param base_vertices: array of lattice vertices, considering it is
        created in a cubic RVE of size 1 and centered on the origin
        :param strut_vertex_pairs: array of strut vertex pairs that define how
        vertices are connected by the struts
        :param cell_size: size of the cubic rve in which the lattice
        cell is enclosed
        :param strut_joints: option to add spherical joints at the vertices
        to better manage strut junctions
        """
        super().__init__(**kwargs)

        self.strut_radius = strut_radius
        self.cell_size = cell_size
        self.strut_joints = strut_joints
        self._strut_heights = strut_heights
        self._base_vertices = base_vertices
        self._strut_vertex_pairs = strut_vertex_pairs

        self.rve = Rve(dim=self.cell_size, center=self.center)

        self.vertices = self._compute_vertices()
        self.strut_centers = self._compute_strut_centers()
        self.strut_directions_cartesian = self._compute_strut_directions()
        self.strut_rotations = self._compute_rotations()

        self._validate_inputs()

        self.cad_shape = self.generate()

    @property
    def base_vertices(self) -> npt.NDArray[np.float64]:
        """Property: coordinates of the vertices for a structure
        centered at the origin and enclosed in a size 1 cubic rve"""
        if self._base_vertices is not None:
            return self._base_vertices
        return self._generate_base_vertices()

    @property
    def strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        """Property: pairs of vertex indices forming a strut"""
        if self._strut_vertex_pairs is not None:
            return self._strut_vertex_pairs
        return self._generate_strut_vertex_pairs()

    @abstractmethod
    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        """Abstract method to generate base vertices, ie as if the
        lattice was centered at the origin and in a cubic size 1 rve.
        """
        pass

    @abstractmethod
    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        """Abstract method to generate strut vertex pairs."""
        pass

    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        return self.center + self.cell_size * self.base_vertices

    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        return np.mean(self.vertices[self.strut_vertex_pairs], axis=1)

    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        vectors = np.diff(self.vertices[self.strut_vertex_pairs], axis=1).squeeze()
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)

    def _validate_inputs(self):
        """Checks coherence of inputs."""

        if self._strut_heights is None:
            raise NotImplementedError("strut_heights must be defined by the subclass")
        if (
            isinstance(self._strut_heights, list)
            and len(self._strut_heights) != self.strut_number
        ):
            raise ValueError(
                f"strut_heights must contain {self.strut_number} values, but {len(self._strut_heights)} were provided."
            )

    @property
    def strut_number(self) -> int:
        return len(self.strut_vertex_pairs)

    @property
    def strut_heights(self) -> list[float]:
        """
        Returns the list of strut lengths.
        If a single value is given, it is converted to a list.
        """
        if isinstance(self._strut_heights, float):
            return [self._strut_heights * self.cell_size] * self.strut_number

        return self._strut_heights * self.cell_size

    def _compute_rotations(self) -> list[Rotation]:
        """Computes rotation from default (1.0, 0.0, 0.0) oriented Cylinder
        for all struts in the lattice using Scipy's Rotation object."""

        default_direction = np.array([1.0, 0.0, 0.0])

        rotations_list = []

        for i in range(self.strut_number):
            if np.all(
                self.strut_directions_cartesian[i] == default_direction
            ) or np.all(self.strut_directions_cartesian[i] == -default_direction):
                rotation_vector = np.zeros(3)
                rotations_list.append(Rotation.from_rotvec(rotation_vector))
            else:
                axis = np.cross(default_direction, self.strut_directions_cartesian[i])
                axis /= np.linalg.norm(axis)
                angle = np.arccos(
                    np.dot(default_direction, self.strut_directions_cartesian[i])
                )
                rotation_vector = angle * axis
                rotations_list.append(Rotation.from_rotvec(rotation_vector))

        return rotations_list

    def generate(self, **_: KwargsGenerateType) -> cq.Shape:
        """Generate a strut-based lattice CAD shape using the given parameters."""
        list_phases: list[Phase] = []
        list_periodic_phases: list[Phase] = []

        for i in range(self.strut_number):
            strut = Cylinder(
                center=tuple(self.strut_centers[i]),
                orientation=self.strut_rotations[i],
                height=self.strut_heights[i],
                radius=self.strut_radius,
            )
            shape = strut.generate()
            list_phases.append(Phase(shape))
        if self.strut_joints:
            for vertex in self.vertices:
                joint = Sphere(
                    center=tuple(vertex),
                    radius=self.strut_radius,
                )
                shape = joint.generate()
                list_phases.append(Phase(shape))

        for phase in list_phases:
            periodic_phase = periodic_split_and_translate(phase=phase, rve=self.rve)
            list_periodic_phases.append(periodic_phase)

        lattice = fuse_shapes(
            [phase.shape for phase in list_periodic_phases],
            retain_edges=False,
        )

        bounding_box = Box(
            center=self.center,
            dim=(self.cell_size, self.cell_size, self.cell_size),
        ).generate()

        cut_lattice = bounding_box.intersect(lattice)

        return cut_lattice

    @property
    def volume(self) -> float:
        volume = self.cad_shape.Volume()

        return volume

    def generate_vtk(
        self,
        size: float = 0.02,
        order: int = 1,
        periodic: bool = True,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate a strut-based lattice VTK shape using the given parameters."""
        cad_lattice = self.cad_shape
        list_phases = [Phase(cad_lattice)]

        with (
            NamedTemporaryFile(suffix=".step", delete=False) as cad_step_file,
            NamedTemporaryFile(suffix=".vtk", delete=False) as mesh_file,
        ):
            cq.exporters.export(cad_lattice, cad_step_file.name)
            if periodic:
                mesh_periodic(
                    mesh_file=cad_step_file.name,
                    rve=self.rve,
                    list_phases=list_phases,
                    size=size,
                    order=order,
                    output_file=mesh_file.name,
                )
            mesh(
                mesh_file=cad_step_file.name,
                list_phases=list_phases,
                size=size,
                order=order,
                output_file=mesh_file.name,
            )

            vtk_lattice = pv.read(mesh_file.name).extract_surface()

        # Solve compatibility issues of NamedTemporaryFiles with Windows
        trash_files_list = [
            cad_step_file.name,
            mesh_file.name,
        ]
        for file in trash_files_list:
            Path(file).unlink()

        return vtk_lattice

    def generateVtk(  # noqa: N802
        self,
        size: float = 0.02,
        order: int = 1,
        periodic: bool = True,
        **kwargs: KwargsGenerateType,
    ) -> pv.PolyData:
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk(
            size=size,
            order=order,
            periodic=periodic,
            **kwargs,
        )
