"""
Phase class to manage list of solids belonging to the same phase
"""
import cadquery as cq


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

        self.centerOfMass = None
        self.inertiaMatrix = None

        Phase.numInstances += 1

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
