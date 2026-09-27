"""
=======================================================================
Abstract Lattice (:mod:`microgen.shape.strut_lattice.abstract_lattice`)
=======================================================================
"""

from __future__ import annotations

from abc import abstractmethod
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.optimize import root_scalar
from scipy.spatial.transform import Rotation

from ...cad import CadShape
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

BALL_POINT_RADIUS_TOLERANCE = 1e-5


class AbstractLattice(Shape):
    """Abstract Class to create strut-based lattice."""

    _UNIT_CUBE_SIZE = 1.0
    _DEFAULT_STRUT_HEIGHTS: float | list[float] | None = None

    def __init__(
        self,
        strut_radius: float | None = None,
        strut_heights: float | list[float] | None = None,
        base_vertices: npt.NDArray[np.float64] | None = None,
        strut_vertex_pairs: npt.NDArray[np.int64] | None = None,
        cell_size: float = 1.0,
        strut_joints: bool = False,
        density: float | None = None,
        **kwargs: Vector3DType | Rotation,
    ) -> None:
        """Abstract Class to create strut-based lattice.

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
        if strut_heights is None:
            strut_heights = type(self)._DEFAULT_STRUT_HEIGHTS
        if strut_radius is not None and density is not None:
            err_msg = (
                "strut radius and density cannot be given at the same time. "
                "Give only one."
            )
            raise ValueError(err_msg)

        if strut_radius is None and density is None:
            err_msg = "strut radius or density must be given. Give one of them."
            raise ValueError(err_msg)

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
        self._cad_shape = None
        self._vtk_shape: tuple[tuple[float, int, bool], pv.PolyData] | None = None

        if density is not None and not 0.0 < density <= 1.0:
            err_msg = f"density must be between 0 and 1. Given: {density}"
            raise ValueError(err_msg)

        self.density = density

        if density is not None:
            self.strut_radius = self._compute_radius_to_fit_density()
        else:
            self.strut_radius = strut_radius

    def _compute_radius_to_fit_density(self) -> float:
        """Solve for the strut radius matching the requested density.

        Each ``root_scalar`` step builds a CAD shape; we stash the final
        one on ``self._cad_shape`` so :meth:`generate_cad` doesn't have to
        rebuild it afterwards.
        """
        RADIUS_MIN = 10e-4
        RADIUS_MAX_MULTIPLIER = 1.0

        def calc_density(radius: float) -> float:
            self.strut_radius = radius
            self._cad_shape = self._generate_cad()
            return self._cad_shape.volume() / (self.cell_size**3)

        return root_scalar(
            lambda radius: float(calc_density(radius)) - self.density,
            bracket=[RADIUS_MIN, RADIUS_MAX_MULTIPLIER * self.cell_size],
        ).root

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
            err_msg = (
                f"strut_heights must contain {self.strut_number} values, "
                f"but {len(self._strut_heights)} were provided."
            )
            raise ValueError(err_msg)

    @property
    def strut_number(self) -> int:
        return len(self.strut_vertex_pairs)

    @property
    def strut_heights(self) -> list[float]:
        """Return the list of strut lengths.

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
                rotation, _ = Rotation.align_vectors(
                    self.strut_directions_cartesian[i], default_direction
                )
                rotations_list.append(rotation)

        return rotations_list

    def generate_cad(self, **_: KwargsGenerateType) -> CadShape:
        if isinstance(self._cad_shape, CadShape):
            return self._cad_shape

        self._cad_shape = self._generate_cad()
        return self._cad_shape

    cad_shape = property(generate_cad)

    def generate_implicit(self: AbstractLattice) -> Shape:
        """Return the lattice as a single composed implicit :class:`Shape`.

        Builds ``(∪ struts ∪ joints) ∩ box`` by F-rep composition — every
        primitive stays as a :class:`Cylinder` / :class:`Sphere` SDF;
        nothing materialises BREP geometry.  Pair with :meth:`Phase.from_shape`
        for an implicit-first pipeline (marching cubes mesh, ``pieces``
        connectivity, grid-quadrature moments — no OCCT required).

        For periodic boundary conditions, the existing CAD path
        (:meth:`generate_cad`) does explicit ``periodic_split_and_translate``;
        an implicit periodic-wrap is not produced here (struts that cross
        the cell boundary are simply clipped by the bounding box).
        """
        from functools import reduce  # noqa: PLC0415
        from operator import or_  # noqa: PLC0415

        primitives: list[Shape] = [
            Cylinder(
                center=tuple(self.strut_centers[i]),
                orientation=self.strut_rotations[i],
                height=self.strut_heights[i],
                radius=self.strut_radius,
            )
            for i in range(self.strut_number)
        ]
        if self.strut_joints:
            primitives.extend(
                Sphere(center=tuple(v), radius=self.strut_radius) for v in self.vertices
            )

        union = reduce(or_, primitives)
        bounding_box = Box(
            center=self.center,
            dim=(self.cell_size, self.cell_size, self.cell_size),
        )
        return union & bounding_box

    def _generate_cad(self, **_: KwargsGenerateType) -> CadShape:
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
            shape = strut.generate_cad()
            list_phases.append(Phase.from_cad(shape))
        if self.strut_joints:
            for vertex in self.vertices:
                joint = Sphere(
                    center=tuple(vertex),
                    radius=self.strut_radius,
                )
                shape = joint.generate_cad()
                list_phases.append(Phase.from_cad(shape))

        for phase in list_phases:
            periodic_phase = periodic_split_and_translate(phase=phase, rve=self.rve)
            list_periodic_phases.append(periodic_phase)

        lattice = fuse_shapes(
            [phase.cad for phase in list_periodic_phases],
            retain_edges=False,
        )

        bounding_box = Box(
            center=self.center,
            dim=(self.cell_size, self.cell_size, self.cell_size),
        ).generate_cad()

        cut_lattice = bounding_box.intersect(lattice)

        return cut_lattice

    @property
    def volume(self) -> float:
        volume = self.cad_shape.volume()

        return volume

    def generate_surface_mesh(
        self,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Return a surface mesh of the lattice (for visualisation).

        Today this delegates to :meth:`mesh_for_fem` with default parameters
        (``size=0.02, order=1, periodic=True``), which runs CAD → STEP →
        gmsh → pyvista. When the F-rep implicit-lattice work lands, this
        method will switch to F-rep marching cubes (no CAD/gmsh required)
        and :meth:`mesh_for_fem` will remain as the explicit FEM-meshing
        path.

        Users who need to control mesh size / element order / periodicity
        should call :meth:`mesh_for_fem` directly.
        """
        return self.mesh_for_fem()

    vtk_shape = property(generate_surface_mesh)

    def mesh_for_fem(
        self,
        size: float = 0.02,
        order: int = 1,
        *,
        periodic: bool = True,
    ) -> pv.PolyData:
        """Build a periodic / non-periodic FEM tet mesh and return its surface.

        Path: ``cad_shape`` → STEP → gmsh (:func:`microgen.mesh_periodic` or
        :func:`microgen.mesh`) → ``pv.read`` → :meth:`extract_surface`.
        Requires the ``[cad]`` extra and gmsh.

        Cached per ``(size, order, periodic)`` tuple on the instance, so
        repeated calls with the same parameters are O(1).

        :param size: target element size (gmsh)
        :param order: element order (gmsh)
        :param periodic: enforce periodicity via :func:`mesh_periodic`
        :return: surface ``pv.PolyData`` extracted from the tet mesh
        """
        params = (size, order, periodic)
        if self._vtk_shape is not None:
            cached_params, cached_mesh = self._vtk_shape
            if cached_params == params:
                return cached_mesh

        cad_lattice = self.cad_shape
        list_phases = [Phase.from_cad(cad_lattice)]

        with (
            NamedTemporaryFile(suffix=".step", delete=False) as cad_step_file,
            NamedTemporaryFile(suffix=".vtk", delete=False) as mesh_file,
        ):
            cad_lattice.export_step(cad_step_file.name)
            mesher = mesh_periodic if periodic else mesh
            mesher_kwargs = (
                {"rve": self.rve, "list_phases": list_phases}
                if periodic
                else {"list_phases": list_phases}
            )
            mesher(
                mesh_file=cad_step_file.name,
                size=size,
                order=order,
                output_file=mesh_file.name,
                **mesher_kwargs,
            )
            vtk_lattice = pv.read(mesh_file.name).extract_surface(algorithm=None)

        # Solve compatibility issues of NamedTemporaryFiles with Windows.
        for tmp in (cad_step_file.name, mesh_file.name):
            Path(tmp).unlink()

        self._vtk_shape = (params, vtk_lattice)
        return vtk_lattice
