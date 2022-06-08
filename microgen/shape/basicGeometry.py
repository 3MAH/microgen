import cadquery as cq

# import box
# from .capsule import Capsule
# from .cylinder import Cylinder
# from .ellipsoid import Ellipsoid
# from .extrudedPolygon import ExtrudedPolygon
# from .polyhedron import Polyhedron
# from .sphere import Sphere
# from .tpms import Tpms


class BasicGeometry:
    """
        BasicGeometry class to manage shapes

        :param shape: name of the shape
        :type shape: str
        :param param_geom: dictionary containing corresponding shape parameters
        :type param_geom: dict
    """

    numInstances = 0

    def __init__(
        self,
        shape: str,
        center: list[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
    ) -> None:

        self.number = self.numInstances
        self.shape = shape
        self.center = center
        self.orientation = orientation
        self.name = self.shape + '_' + str(self.number)

        self.geometry = None  # type: cq.Shape
        BasicGeometry.numInstances += 1
