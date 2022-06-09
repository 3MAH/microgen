"""
Phase class to manage list of solids belonging to the same phase
"""

import numpy as np
import cadquery as cq
from OCP.gp import (gp_Vec,gp_Mat)
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp

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

        self.shape = shape
        self.solids = solids
        self.center = center
        self.orientation = orientation

        self.name = "Phase_" + str(self.numInstances)

        self._centerOfMass = None
        self._inertiaMatrix = None

        Phase.numInstances += 1

    @property
    def centerOfMass(self):
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
        BRepGProp.VolumeProperties_s(self.shape.wrapped, Properties)
 
        COM = Properties.CentreOfMass()
        self._centerOfMass = np.array([COM.X(),COM.Y(),COM.Z()])
        
    @property
    def inertiaMatrix(self):
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
        BRepGProp.VolumeProperties_s(self.shape.wrapped, Properties)

        INM = Properties.MatrixOfInertia()
        self._inertiaMatrix = np.array([[INM.Value(1,1), INM.Value(1,2), INM.Value(1,3)], [INM.Value(2,1), INM.Value(2,2), INM.Value(2,3)], [INM.Value(3,1), INM.Value(3,2), INM.Value(3,3)]])

    def getSolids(self) -> list[cq.Solid]:
        if len(self.solids) > 0:
            return self.solids
        elif self.shape is not None:
            return self.shape.Solids()
        else:
            print("No solids")

    def getFlatSolidList(self) -> list[cq.Solid]:
        if isinstance(self.solids[0], list):  # if solids is list of list
            return [item.copy() for sublist in self.solids for item in sublist]
        else:
            return self.solids
