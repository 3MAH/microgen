"""Phase class to manage list of solids belonging to the same phase."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Sequence

import cadquery as cq
import numpy as np
import numpy.typing as npt
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps

if TYPE_CHECKING:
    from .rve import Rve


class Phase:
    """Phase class to manage list of solids belonging to the same phase properties.

    - centerOfMass
    - inertiaMatrix
    - shape
    - solids

    :param shape: Shape object
    :param solids: list of cq.Solid or list of list
    :param center: center
    :param orientation: orientation
    """

    num_instances = 0

    def __init__(
        self: Phase,
        shape: cq.Shape | None = None,
        solids: list[cq.Solid] | None = None,
        center: tuple[float, float, float] | None = None,
        orientation: tuple[float, float, float] | None = None,
    ) -> None:
        """Initialize the phase object."""
        self._shape = shape
        self._solids: list[cq.Solid] = solids if solids is not None else []
        self.center = center
        self.orientation = orientation

        if shape is None and solids == []:
            warnings.warn("Empty phase", stacklevel=2)

        self.name = f"Phase_{self.num_instances}"

        self._center_of_mass = None
        self._inertia_matrix = None

        Phase.num_instances += 1

    def get_center_of_mass(
        self: Phase,
        *,
        compute: bool = True,
    ) -> npt.NDArray[np.float64]:
        """Return the center of 'mass' of an object.

        :param compute: if False and centerOfMass already exists, \
            does not compute it (use carefully)
        """
        if isinstance(self._center_of_mass, np.ndarray) and not compute:
            return self._center_of_mass

        self._compute_center_of_mass()
        return self._center_of_mass

    center_of_mass = property(get_center_of_mass)

    def _compute_center_of_mass(self: Phase) -> None:
        """Calculate the center of 'mass' of an object."""
        properties = GProp_GProps()
        BRepGProp.VolumeProperties_s(self._shape.wrapped, properties)

        com = properties.CentreOfMass()
        self._center_of_mass = np.array([com.X(), com.Y(), com.Z()])

    def get_inertia_matrix(
        self: Phase,
        *,
        compute: bool = True,
    ) -> npt.NDArray[np.float64]:
        """Calculate the inertia Matrix of an object.

        :param compute: if False and inertiaMatrix already exists, \
            does not compute it (use carefully)
        """
        if isinstance(self._inertia_matrix, np.ndarray) and not compute:
            return self._inertia_matrix

        self._compute_inertia_matrix()
        return self._inertia_matrix

    inertia_matrix = property(get_inertia_matrix)

    def _compute_inertia_matrix(self: Phase) -> None:
        """Calculate the inertia Matrix of an object."""
        properties = GProp_GProps()
        BRepGProp.VolumeProperties_s(self._shape.wrapped, properties)

        inm = properties.MatrixOfInertia()
        self._inertia_matrix = np.array(
            [
                [inm.Value(1, 1), inm.Value(1, 2), inm.Value(1, 3)],
                [inm.Value(2, 1), inm.Value(2, 2), inm.Value(2, 3)],
                [inm.Value(3, 1), inm.Value(3, 2), inm.Value(3, 3)],
            ],
        )

    @property
    def shape(self: Phase) -> cq.Shape | None:
        """Return the shape of the phase."""
        if self._shape is not None:
            return self._shape
        if len(self._solids) > 0:
            # there may be a faster way
            compound = cq.Compound.makeCompound(self._solids)
            self._shape = cq.Shape(compound.wrapped)
            return self._shape

        warnings.warn("No shape or solids", stacklevel=2)
        return None

    @property
    def solids(self: Phase) -> list[cq.Solid]:
        """Return the list of solids of the phase."""
        if len(self._solids) > 0:
            return self._solids
        if self._shape is not None:
            self._solids = self._shape.Solids()
            return self._solids

        warnings.warn("No solids or shape", stacklevel=2)
        return []

    def translate(self: Phase, vec: Sequence[float]) -> None:
        """Translate phase by a given vector."""
        self._shape.move(cq.Location(cq.Vector(*vec)))
        self._compute_center_of_mass()

    @staticmethod
    def rescale_shape(
        shape: cq.Shape,
        scale: float | tuple[float, float, float],
    ) -> cq.Shape:
        """Rescale given object according to scale parameters [dim_x, dim_y, dim_z].

        :param shape: Shape
        :param scale: float or list of scale factor in each direction

        :return shape: rescaled Shape
        """
        if isinstance(scale, float):
            scale = (scale, scale, scale)

        center = shape.Center()

        # move the shape at (0, 0, 0) to rescale it
        shape.move(cq.Location(cq.Vector(-center.x, -center.y, -center.z)))

        # then move it back to its center with transform Matrix
        transform_mat = cq.Matrix(
            [
                [scale[0], 0, 0, center.x],
                [0, scale[1], 0, center.y],
                [0, 0, scale[2], center.z],
            ],
        )

        return shape.transformGeometry(transform_mat)

    def rescale(self: Phase, scale: float | tuple[float, float, float]) -> None:
        """Rescale phase according to scale parameters [dim_x, dim_y, dim_z].

        :param scale: float or list of scale factor in each direction
        """
        self._shape = self.rescale_shape(self._shape, scale)

    @staticmethod
    def repeat_shape(
        unit_geom: cq.Shape,
        rve: Rve,
        grid: tuple[int, int, int],
    ) -> cq.Shape:
        """Repeat unit geometry in each direction according to the given grid.

        :param unit_geom: Shape to repeat
        :param rve: RVE of the geometry to repeat
        :param grid: list of number of geometry repetitions in each direction

        :return: cq shape of the repeated geometry
        """
        center = np.array(unit_geom.Center().toTuple())

        xyz_repeat = cq.Assembly()
        for idx in np.ndindex(*grid):
            pos = center - rve.dim * (0.5 * np.array(grid) - 0.5 - np.array(idx))
            xyz_repeat.add(unit_geom, loc=cq.Location(cq.Vector(*pos)))

        return cq.Shape(xyz_repeat.toCompound().wrapped)

    def repeat(self: Phase, rve: Rve, grid: tuple[int, int, int]) -> None:
        """Repeat phase in each direction according to the given grid.

        :param rve: RVE of the phase to repeat
        :param grid: list of number of phase repetitions in each direction
        """
        self._shape = self.repeat_shape(self.shape, rve, grid)

    def split_solids(self: Phase, rve: Rve, grid: list[int]) -> list[cq.Solid]:
        """Split solids from phase according to the rve divided by the given grid.

        :param rve: RVE divided by the given grid
        :param grid: number of divisions in each direction [x, y, z]

        :return: list of solids
        """
        solids: list[cq.Solid] = []

        for solid in self.solids:
            wk_plane = cq.Workplane().add(solid)
            for dim in range(3):
                direction = cq.Vector([int(dim == i) for i in range(3)])
                coords = np.linspace(
                    start=rve.min_point[dim],
                    stop=rve.max_point[dim],
                    num=grid[dim],
                    endpoint=False,
                )[1:]
                for pos in coords:
                    point = pos * direction
                    plane = cq.Face.makePlane(basePnt=point, dir=direction)
                    wk_plane = wk_plane.split(plane)

            solids += wk_plane.val().Solids()
        return solids

    def rasterize(
        self: Phase,
        rve: Rve,
        grid: list[int],
        phasePerRaster: bool | None = None,  # noqa: N803
        *,
        phase_per_raster: bool = True,
    ) -> list[Phase] | None:
        """Raster solids from phase according to the rve divided by the given grid.

        :param rve: RVE divided by the given grid
        :param grid: number of divisions in each direction [x, y, z]
        :param phase_per_raster: if True, returns list of phases

        :return: list of Phases if required
        """
        if phasePerRaster is not None:
            warnings.warn(
                "phasePerRaster is deprecated, use phase_per_raster instead",
                DeprecationWarning,
                stacklevel=2,
            )
            phase_per_raster = phasePerRaster

        solids: list[cq.Solid] = self.split_solids(rve, grid)

        if phase_per_raster:
            return self.generate_phase_per_raster(solids, rve, grid)
        self._solids = solids
        compound = cq.Compound.makeCompound(self._solids)
        self._shape = cq.Shape(compound.wrapped)
        return None

    @classmethod
    def generate_phase_per_raster(
        cls: type[Phase],
        solids: list[cq.Solid],
        rve: Rve,
        grid: list[int],
    ) -> list[Phase]:
        """Raster solids from phase according to the rve divided by the given grid.

        :param solidList: list of solids
        :param grid: number of divisions in each direction [x, y, z]
        :param rve: RVE divided by the given grid

        :return: list of Phases
        """
        grid = np.array(grid)
        solids_phases: list[list[cq.Solid]] = [[] for _ in range(np.prod(grid))]
        for solid in solids:
            center = np.array(solid.Center().toTuple())
            i, j, k = np.floor(grid * (center - rve.min_point) / rve.dim).astype(int)
            ind = i + grid[0] * j + grid[0] * grid[1] * k
            solids_phases[ind].append(solid)
        return [Phase(solids=solids) for solids in solids_phases if len(solids) > 0]

    # Deprecated methods

    @classmethod
    def generatePhasePerRaster(  # noqa: N802
        cls: type[Phase],
        solidList: list[cq.Solid],  # noqa: N803
        rve: Rve,
        grid: list[int],
    ) -> list[Phase]:
        """See generate_phase_per_raster.

        Deprecated in favor of generate_phase_per_raster.
        """
        warnings.warn(
            "generatePhasePerRaster is deprecated, \
                use generate_phase_per_raster instead",
            DeprecationWarning,
            stacklevel=2,
        )
        return cls.generate_phase_per_raster(solidList, rve, grid)

    @staticmethod
    def repeatShape(  # noqa: N802
        unit_geom: cq.Shape,
        rve: Rve,
        grid: tuple[int, int, int],
    ) -> cq.Shape:
        """See repeat_shape.

        Deprecated in favor of repeat_shape.
        """
        warnings.warn(
            "repeatShape is deprecated, use repeat_shape instead",
            DeprecationWarning,
            stacklevel=2,
        )
        return Phase.repeat_shape(unit_geom, rve, grid)

    @staticmethod
    def rescaleShape(  # noqa: N802
        shape: cq.Shape,
        scale: float | tuple[float, float, float],
    ) -> cq.Shape:
        """See rescale_shape.

        Deprecated in favor of rescale_shape.
        """
        warnings.warn(
            "rescaleShape is deprecated, use rescale_shape instead",
            DeprecationWarning,
            stacklevel=2,
        )
        return Phase.rescale_shape(shape, scale)

    def getCenterOfMass(self: Phase, compute: bool = True) -> np.ndarray:  # noqa: N802, FBT001, FBT002
        """See get_center_of_mass.

        Deprecated in favor of get_center_of_mass.
        """
        warnings.warn(
            "centerOfMass is deprecated, use center_of_mass instead",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.get_center_of_mass(compute=compute)

    centerOfMass = property(getCenterOfMass)  # noqa: N815

    def getInertiaMatrix(self: Phase, compute: bool = True) -> np.ndarray:  # noqa: N802, FBT001, FBT002
        """See get_inertia_matrix.

        Deprecated in favor of get_inertia_matrix.
        """
        warnings.warn(
            "inertiaMatrix is deprecated, use inertia_matrix instead",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.get_inertia_matrix(compute=compute)

    inertiaMatrix = property(getInertiaMatrix)  # noqa: N815
