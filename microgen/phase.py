"""
Phase class to manage list of solids belonging to the same phase
"""

import cadquery as cq
import numpy as np
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps

from typing import Union


class Phase:
    """
        Phase class to manage list of solids belonging to the same phase

    :param shape: Shape object
    :param solids: list of cq.Solid or list of list
    :param center: center
    :param orientation: orientation
    """

    numInstances = 0

    def __init__(
        self,
        shape: cq.Shape = None,
        solids: list[cq.Solid] = [],
        center: tuple[float, float, float] = None,
        orientation: tuple[float, float, float] = None,
    ) -> None:

        if shape is None and solids == []:
            print("Empty phase")

        self._shape = shape
        self._solids = solids
        self.center = center
        self.orientation = orientation

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

    def _computeCenterOfMass(self):
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

    def _computeInertiaMatrix(self):
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
    def shape(self) -> Union[cq.Shape, None]:
        if self._shape is not None:
            return self._shape
        elif len(self._solids) > 0:
            # there may be a fastest way
            compound = cq.Compound.makeCompound(self._solids)
            self._shape = cq.Shape(compound.wrapped)
            return self._shape
        else:
            print("No shape or solids")
            return None

    @property
    def solids(self) -> list[cq.Solid]:
        if len(self._solids) > 0:
            return self._solids
        elif self._shape is not None:
            self._solids = self._shape.Solids()
            return self._solids
        else:
            print("No solids or shape")
            return []

    def translate(self, vec: Union[tuple, np.ndarray]) -> None:
        if type(vec) == tuple:
            self._shape.move(cq.Location(cq.Vector(vec)))
        else:
            self._shape.move(cq.Location(cq.Vector(vec[0], vec[1], vec[2])))
        self._computeCenterOfMass()

    centerOfMass = property(getCenterOfMass)
    inertiaMatrix = property(getInertiaMatrix)
