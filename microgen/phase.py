"""
Phase class to manage list of solids belonging to the same phase
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import cadquery as cq
import numpy as np
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps

from .rve import Rve


class Phase:
    """
    Phase class to manage list of solids belonging to the same phase
    properties: centerOfMass, inertiaMatrix, shape, solids

    :param shape: Shape object
    :param solids: list of cq.Solid or list of list
    :param center: center
    :param orientation: orientation
    """

    numInstances = 0

    def __init__(
        self,
        shape: Optional[cq.Shape] = None,
        solids: Optional[List[cq.Solid]] = None,
        center: Optional[Tuple[float, float, float]] = None,
        orientation: Optional[Tuple[float, float, float]] = None,
    ) -> None:
        self._shape = shape
        self._solids: List[cq.Solid] = solids if solids is not None else []
        self.center = center
        self.orientation = orientation

        if shape is None and solids == []:
            print("Empty phase")

        self.name = "Phase_" + str(self.numInstances)

        self._centerOfMass = None
        self._inertiaMatrix = None

        Phase.numInstances += 1

    def getCenterOfMass(self, compute: bool = True) -> np.ndarray:
        """
        Returns the center of 'mass' of an object.
        :param compute: if False and centerOfMass already exists, does not compute it (use carefully)
        """
        if isinstance(self._centerOfMass, np.ndarray) and not compute:
            return self._centerOfMass
        else:
            self._computeCenterOfMass()
            return self._centerOfMass

    centerOfMass = property(getCenterOfMass)

    def _computeCenterOfMass(self) -> None:
        """
        Calculates the center of 'mass' of an object.
        """
        Properties = GProp_GProps()
        BRepGProp.VolumeProperties_s(self._shape.wrapped, Properties)

        com = Properties.CentreOfMass()
        self._centerOfMass = np.array([com.X(), com.Y(), com.Z()])

    def getInertiaMatrix(self, compute: bool = True) -> np.ndarray:
        """
        Calculates the inertia Matrix of an object.
        :param compute: if False and inertiaMatrix already exists, does not compute it (use carefully)
        """
        if isinstance(self._inertiaMatrix, np.ndarray) and not compute:
            return self._inertiaMatrix
        else:
            self._computeInertiaMatrix()
            return self._inertiaMatrix

    inertiaMatrix = property(getInertiaMatrix)

    def _computeInertiaMatrix(self) -> None:
        """
        Calculates the inertia Matrix of an object.
        """
        Properties = GProp_GProps()
        BRepGProp.VolumeProperties_s(self._shape.wrapped, Properties)

        inm = Properties.MatrixOfInertia()
        self._inertiaMatrix = np.array(
            [
                [inm.Value(1, 1), inm.Value(1, 2), inm.Value(1, 3)],
                [inm.Value(2, 1), inm.Value(2, 2), inm.Value(2, 3)],
                [inm.Value(3, 1), inm.Value(3, 2), inm.Value(3, 3)],
            ]
        )

    @property
    def shape(self) -> Optional[cq.Shape]:
        if self._shape is not None:
            return self._shape
        elif len(self._solids) > 0:
            # there may be a faster way
            compound = cq.Compound.makeCompound(self._solids)
            self._shape = cq.Shape(compound.wrapped)
            return self._shape
        else:
            print("No shape or solids")
            return None

    @property
    def solids(self) -> List[cq.Solid]:
        if len(self._solids) > 0:
            return self._solids
        elif self._shape is not None:
            self._solids = self._shape.Solids()
            return self._solids
        else:
            print("No solids or shape")
            return []

    def translate(self, vec: Sequence[float]) -> None:
        self._shape.move(cq.Location(cq.Vector(*vec)))
        self._computeCenterOfMass()

    @staticmethod
    def rescaleShape(
        shape: cq.Shape, scale: float | Tuple[float, float, float]
    ) -> cq.Shape:
        """
        Rescale given object according to scale parameters [dim_x, dim_y, dim_z]

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
            ]
        )
        shape = shape.transformGeometry(transform_mat)

        return shape

    def rescale(self, scale: float | Tuple[float, float, float]) -> None:
        """
        Rescale phase according to scale parameters [dim_x, dim_y, dim_z]

        :param scale: float or list of scale factor in each direction
        """
        self._shape = self.rescaleShape(self._shape, scale)

    @staticmethod
    def repeatShape(
        unit_geom: cq.Shape, rve: Rve, grid: Tuple[int, int, int]
    ) -> cq.Shape:
        """
        Repeats unit geometry in each direction according to the given grid

        :param unit_geom: Shape to repeat
        :param rve: RVE of the geometry to repeat
        :param grid: list of number of geometry repetitions in each direction

        :return: cq shape of the repeated geometry
        """

        center = unit_geom.Center()

        xyz_repeat = cq.Assembly()
        for i_x in range(grid[0]):
            for i_y in range(grid[1]):
                for i_z in range(grid[2]):
                    xyz_repeat.add(
                        unit_geom,
                        loc=cq.Location(
                            cq.Vector(
                                center.x - rve.dim_x * (0.5 * grid[0] - 0.5 - i_x),
                                center.y - rve.dim_y * (0.5 * grid[1] - 0.5 - i_y),
                                center.z - rve.dim_z * (0.5 * grid[2] - 0.5 - i_z),
                            )
                        ),
                    )
        compound = xyz_repeat.toCompound()
        shape = cq.Shape(compound.wrapped)
        return shape

    def repeat(self, rve: Rve, grid: Tuple[int, int, int]) -> None:
        """
        Repeats phase in each direction according to the given grid

        :param rve: RVE of the phase to repeat
        :param grid: list of number of phase repetitions in each direction
        """
        self._shape = self.repeatShape(self.shape, rve, grid)

    def buildSolids(self, rve: Rve, grid: List[int]) -> List[cq.Solid]:
        """
        Build solids from phase according to the rve divided by the given grid

        :param rve: RVE divided by the given grid
        :param grid: number of divisions in each direction [x, y, z]

        :return: list of solids
        """
        solidList: List[cq.Solid] = []

        for solid in self.solids:
            wk_plane = cq.Workplane().add(solid)
            xgrid = np.linspace(rve.x_min, rve.x_max, num=grid[0])
            ygrid = np.linspace(rve.y_min, rve.y_max, num=grid[1])
            zgrid = np.linspace(rve.z_min, rve.z_max, num=grid[2])
            np.delete(xgrid, 0)
            np.delete(ygrid, 0)
            np.delete(zgrid, 0)
            for i in xgrid:
                Plane_x = cq.Face.makePlane(basePnt=(i, 0, 0), dir=(1, 0, 0))
                wk_plane = wk_plane.split(cq.Workplane().add(Plane_x))
            for j in ygrid:
                Plane_y = cq.Face.makePlane(basePnt=(0, j, 0), dir=(0, 1, 0))
                wk_plane = wk_plane.split(cq.Workplane().add(Plane_y))
            for k in zgrid:
                Plane_z = cq.Face.makePlane(basePnt=(0, 0, k), dir=(0, 0, 1))
                wk_plane = wk_plane.split(cq.Workplane().add(Plane_z))

            for subsolid in wk_plane.val().Solids():
                solidList.append(subsolid)
        return solidList

    def rasterize(
        self, rve: Rve, grid: List[int], phasePerRaster: bool = True
    ) -> Optional[List[Phase]]:
        """
        Rasters solids from phase according to the rve divided by the given grid

        :param rve: RVE divided by the given grid
        :param grid: number of divisions in each direction [x, y, z]
        :param phasePerRaster: if True, returns list of phases

        :return: list of Phases if required
        """
        solidList: List[cq.Solid] = self.buildSolids(rve, grid)

        if phasePerRaster:
            return self.generatePhasePerRaster(solidList, rve, grid)
        self._solids = solidList
        compound = cq.Compound.makeCompound(self._solids)
        self._shape = cq.Shape(compound.wrapped)
        return None

    @classmethod
    def generatePhasePerRaster(
        cls, solidList: List[cq.Solid], rve: Rve, grid: List[int]
    ) -> List[Phase]:
        """
        Rasters solids from phase according to the rve divided by the given grid

        :param solidList: list of solids
        :param grid: number of divisions in each direction [x, y, z]
        :param rve: RVE divided by the given grid

        :return: list of Phases
        """
        solids_phases: List[List[cq.Solid]] = [
            [] for _ in range(grid[0] * grid[1] * grid[2])
        ]
        for solid in solidList:
            center = solid.Center()
            i = int(round((center.x - rve.x_min) / (rve.dx / grid[0])))
            j = int(round((center.y - rve.y_min) / (rve.dy / grid[1])))
            k = int(round((center.z - rve.z_min) / (rve.dz / grid[2])))
            ind = i + grid[0] * j + grid[0] * grid[1] * k
            solids_phases[ind].append(solid)
        return [cls(solids=solids) for solids in solids_phases if len(solids) > 0]
