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
        solids: list = [],
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

    @property
    def centerOfMass(self) -> np.ndarray:
        """
        Returns the center of 'mass' of an object.
        """
        if isinstance(self._centerOfMass, np.ndarray):
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

    @property
    def inertiaMatrix(self) -> np.ndarray:
        """
        Calculates the inertia Matrix of an object.
        """
        if isinstance(self._inertiaMatrix, np.ndarray):
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

    def getFlatSolidList(self) -> list[cq.Solid]:
        if isinstance(self._solids[0], list):  # if solids is list of list
            return [item.copy() for sublist in self._solids for item in sublist]
        else:
            return self._solids
