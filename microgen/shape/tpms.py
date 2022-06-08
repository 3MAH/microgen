"""
=============================================
TPMS (:mod:`microgen.shape.tpms`)
=============================================
"""
import math
import os
from typing import Callable

import cadquery as cq
import pygalmesh
from numpy import cos, pi, sin
from OCP.StlAPI import StlAPI_Reader
from OCP.TopoDS import TopoDS_Shape

from .basicGeometry import BasicGeometry
from ..operations import fuseParts


class Tpms(BasicGeometry):
    """
    Class to generate Triply Periodical Minimal Surfaces (TPMS)
    geometry from a given mathematical function, with given thickness

    functions available :
        - :class:`~microgen.shape.tpms.gyroid`
        - :class:`~microgen.shape.tpms.schwarzP`
        - :class:`~microgen.shape.tpms.schwarzD`
        - :class:`~microgen.shape.tpms.neovius`
        - :class:`~microgen.shape.tpms.schoenIWP`
        - :class:`~microgen.shape.tpms.schoenFRD`
        - :class:`~microgen.shape.tpms.fischerKochS`
        - :class:`~microgen.shape.tpms.pmy`
        - :class:`~microgen.shape.tpms.honeycomb`
    """
    def __init__(
        self,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        surface_function: Callable[[float, float, float, float], float] = None,
        type_part: str = 'sheet',
        thickness: float = 0,
        path_data: str = ".",
    ) -> None:
        super().__init__(shape='TPMS', center=center, orientation=orientation)
        self.surface_function = surface_function
        self.type_part = type_part
        self.thickness = thickness
        self.path_data = path_data

    def createSurface(
        self,
        isovalue: float = 0,
        filename: str = "surface.stl",
        sizeMesh: float = 0.05,
        minFacetAngle: float = 10,
        maxRadius: float = 0.05,
        verbose: bool = False,
    ) -> None:
        surface = Generator(isovalue, self.surface_function)
        mesh = pygalmesh.generate_surface_mesh(
            surface,
            min_facet_angle=minFacetAngle,
            max_radius_surface_delaunay_ball=maxRadius,
            max_facet_distance=sizeMesh,
            verbose=verbose,
        )
        if not os.path.isdir(self.path_data):
            os.mkdir(self.path_data)
        mesh.write(self.path_data + "/" + filename)

    def generateSurface(
        self,
        sizeMesh: float = 0.05,
        minFacetAngle: float = 10,
        maxRadius: float = 0.05,
        verbose: bool = False,
    ) -> cq.Shape:
        """
        Creates TPMS sheet surface
        """
        self.createSurface(0, "surface.stl", sizeMesh, minFacetAngle, maxRadius, verbose)
        ocp_shape = TopoDS_Shape()
        stl_reader = StlAPI_Reader()
        stl_reader.Read(ocp_shape, self.path_data + "/" + "surface.stl")

        box = cq.Shape(cq.Workplane("front").box(1, 1, 1).val().wrapped)
        surface = cq.Shape(ocp_shape)
        return surface.intersect(box)

    def generate(
        self,
        sizeMesh: float = 0.05,
        minFacetAngle: float = 10,
        maxRadius: float = 0.05,
        verbose: bool = False,
    ) -> cq.Shape:
        """
        Creates thick TPMS geometry (sheet or skeletal part) from surface
        files generated in the directory given by path_data
        """

        thickness = self.thickness * pi
        self.createSurface(thickness / 4.0, "tpms_testplus.stl", sizeMesh, minFacetAngle, maxRadius, verbose)
        self.createSurface(-thickness / 4.0, "tpms_testminus.stl", sizeMesh, minFacetAngle, maxRadius, verbose)
        self.createSurface(thickness / 2.0, "tpms_plus.stl", sizeMesh, minFacetAngle, maxRadius, verbose)
        self.createSurface(-thickness / 2.0, "tpms_minus.stl", sizeMesh, minFacetAngle, maxRadius, verbose)

        surf_tp = TopoDS_Shape()
        surf_tm = TopoDS_Shape()
        surf_p = TopoDS_Shape()
        surf_m = TopoDS_Shape()
        stl_reader = StlAPI_Reader()

        if not (
            stl_reader.Read(surf_tp, self.path_data + "/" + "tpms_testplus.stl")
            and stl_reader.Read(surf_tm, self.path_data + "/" + "tpms_testminus.stl")
            and stl_reader.Read(surf_p, self.path_data + "/" + "tpms_plus.stl")
            and stl_reader.Read(surf_m, self.path_data + "/" + "tpms_minus.stl")
        ):
            raise ValueError(
                "tpms_plus.stl, tpms_minus.stl, tpms_testplus.stl, tpms_testminus.stl files not found in '"
                + self.path_data
                + "' folder"
            )

        face_cut_tp = cq.Face(surf_tp)
        face_cut_tm = cq.Face(surf_tm)
        face_cut_p = cq.Face(surf_p)
        face_cut_m = cq.Face(surf_m)

        box = cq.Workplane("front").box(1, 1, 1)

        boxCut = box.split(face_cut_p)
        boxCut = boxCut.split(face_cut_m)

        boxSolids = boxCut.solids().all()  # type: list[cq.Workplane]

        listSolids = []  # type: list[tuple[int, cq.Shape]]

        for solid in boxSolids:
            temp = solid.split(face_cut_tp)
            temp = temp.split(face_cut_tm)
            tempSize = temp.solids().size()
            listSolids.append((tempSize, solid.val()))

        sheet = [el[1] for el in listSolids if el[0] > 1]
        skeletal = [el[1] for el in listSolids if el[0] == 1]

        if self.type_part == "sheet":
            to_fuse = [cq.Shape(s.wrapped) for s in sheet]
            return_object = fuseParts(to_fuse, True)
        elif self.type_part == "skeletal":
            to_fuse = [cq.Shape(s.wrapped) for s in skeletal]
            return_object = fuseParts(to_fuse, False)
        return return_object.shape


class Generator(pygalmesh.DomainBase):
    """
    Base class to generate TPMS geometry with a given surface function
    """
    def __init__(
        self,
        height: float,
        surface_function: Callable[[float, float, float, float], float],
    ) -> None:
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = 1
        self.waist_radius = math.sqrt(0.5 ** 2 + 0.5 ** 2)
        self.bounding_sphere_squared_radius = (
            math.sqrt(0.5 ** 2 + 0.5 ** 2 + 0.5 ** 2)
            * 1.1
        )
        self.surface_function = surface_function

    def get_bounding_sphere_squared_radius(self) -> float:
        return self.bounding_sphere_squared_radius

    def eval(self, pos: list) -> float:
        return self.surface_function(pos[0], pos[1], pos[2], self.height)


# Lidinoid -> 0.5*(sin(2*x)*cos(y)*sin(z) + sin(2*y)*cos(z)*sin(x) + sin(2*z)*cos(x)*sin(y)) - 0.5*(cos(2*x)*cos(2*y) + cos(2*y)*cos(2*z) + cos(2*z)*cos(2*x)) + 0.15 = 0


def gyroid(x: float, y: float, z: float, height: float) -> float:
    """
    :math:`sin(2 \pi x) cos(2 \pi y) + sin(2 \pi y) cos(2 \pi z) + sin(2 \pi z) cos(2 \pi x) = 0`
    """
    if abs(x) + abs(y) + abs(z) > 1.0e-8:
        return (
            sin(2 * pi * x) * cos(2 * pi * y)
            + sin(2 * pi * y) * cos(2 * pi * z)
            + sin(2 * pi * z) * cos(2 * pi * x)
            + height
        )
    else:
        return 1.0


def schwarzP(x: float, y: float, z: float, height: float) -> float:
    """
    :math:`cos(2 \pi x) + cos(2 \pi y) + cos(2 \pi z) = 0`
    """
    return cos(2 * pi * x) + cos(2 * pi * y) + cos(2 * pi * z) + height


def schwarzD(x: float, y: float, z: float, height: float) -> float:
    """
    :math:`sin(2 \pi x) sin(2 \pi y) sin(2 \pi z) + \
           sin(2 \pi x) cos(2 \pi y) cos(2 \pi z) + \
           cos(2 \pi x) sin(2 \pi y) cos(2 \pi z) + \
           cos(2 \pi x) cos(2 \pi y) sin(2 \pi z) = 0`
    """
    a = sin(2 * pi * x) * sin(2 * pi * y) * sin(2 * pi * z)
    b = sin(2 * pi * x) * cos(2 * pi * y) * cos(2 * pi * z)
    c = cos(2 * pi * x) * sin(2 * pi * y) * cos(2 * pi * z)
    d = cos(2 * pi * x) * cos(2 * pi * y) * sin(2 * pi * z)
    return a + b + c + d + height


def neovius(x: float, y: float, z: float, height: float) -> float:
    """
    :math:`3 cos(2 \pi x) + cos(2 \pi y) + cos(2 \pi z) + \
           4 cos(2 \pi x) cos(2 \pi y) cos(2 \pi z) = 0`
    """
    a = 3 * cos(2 * pi * x) + cos(2 * pi * y) + cos(2 * pi * z)
    b = 4 * cos(2 * pi * x) * cos(2 * pi * y) * cos(2 * pi * z)

    return a + b + height


def schoenIWP(x: float, y: float, z: float, height: float) -> float:
    """
    :math:`2 ( cos(2 \pi x) cos(2 \pi y) + \
               cos(2 \pi y) cos(2 \pi z) + \
               cos(2 \pi z) cos(2 \pi x)) - \
           (cos(4 \pi x) + cos(4 \pi y) + cos(4 \pi z)) = 0`
    """
    a = 2 * (
        cos(2 * pi * x) * cos(2 * pi * y)
        + cos(2 * pi * y) * cos(2 * pi * z)
        + cos(2 * pi * z) * cos(2 * pi * x)
    )
    b = cos(4 * pi * x) + cos(4 * pi * y) + cos(4 * pi * z)

    return a - b + height


def schoenFRD(x: float, y: float, z: float, height: float) -> float:
    """
    :math:`4 cos(2 \pi x) cos(2 \pi y) cos(2 \pi z) - \
           (cos(4 \pi x) cos(4 \pi y) + \
            cos(4 \pi y) cos(4 \pi z) + \
            cos(4 \pi z) cos(4 \pi x)) = 0`
    """
    a = 4 * cos(2 * pi * x) * cos(2 * pi * y) * cos(2 * pi * z)
    b = (
        cos(4 * pi * x) * cos(4 * pi * y)
        + cos(4 * pi * y) * cos(4 * pi * z)
        + cos(4 * pi * z) * cos(4 * pi * x)
    )
    return a - b + height


def fischerKochS(x: float, y: float, z: float, height: float) -> float:
    """
    :math:`cos(4 \pi x) sin(2 \pi y) cos(2 \pi z) + \
           cos(2 \pi x) cos(4 \pi y) sin(2 \pi z) + \
           sin(2 \pi x) cos(2 \pi y) cos(4 \pi z) = 0`
    """
    a = cos(4 * pi * x) * sin(2 * pi * y) * cos(2 * pi * z)
    b = cos(2 * pi * x) * cos(4 * pi * y) * sin(2 * pi * z)
    c = sin(2 * pi * x) * cos(2 * pi * y) * cos(4 * pi * z)

    return a + b + c + height


def pmy(x: float, y: float, z: float, height: float) -> float:
    """
    :math:`2 cos(2 \pi x) cos(2 \pi y) cos(2 \pi z) + \
           sin(4 \pi x) sin(2 \pi y) + \
           sin(2 \pi x) sin(4 \pi z) + \
           sin(4 \pi y) sin(2 \pi z) = 0`
    """
    a = 2 * cos(2 * pi * x) * cos(2 * pi * y) * cos(2 * pi * z)
    b = sin(4 * pi * x) * sin(2 * pi * y)
    c = sin(2 * pi * x) * sin(4 * pi * z)
    d = sin(4 * pi * y) * sin(2 * pi * z)

    return a + b + c + d + height


def honeycomb(x: float, y: float, z: float, height: float) -> float:
    """
    :math:`sin(2 \pi x) cos(2 \pi y) + sin(2 \pi y) + cos(2 \pi z) = 0`
    """
    return (
        sin(2 * pi * x) * cos(2 * pi * y) + sin(2 * pi * y) + cos(2 * pi * z) + height
    )
