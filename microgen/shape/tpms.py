"""
=============================================
TPMS (:mod:`microgen.shape.tpms`)
=============================================
"""
import math
import os
from typing import Callable

import cadquery as cq
import numpy as np
import pygalmesh
from numpy import cos, pi, sin
from OCP.StlAPI import StlAPI_Reader
from OCP.TopoDS import TopoDS_Shape

from ..operations import fuseParts
from ..rve import Rve


class Tpms:
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
        center: np.ndarray,
        orientation: np.ndarray,
        surface_function: Callable[[float, float, float, float], float],
        type_part: str,
        thickness: float,
        number: int = 0,
    ) -> None:
        self.center = center
        self.orientation = orientation
        self.surface_function = surface_function
        self.type_part = type_part
        self.thickness = thickness
        self.number = number

        self.name_part = "tpms" + str(self.number)

    def createSurfaces(
        self,
        rve: Rve,
        sizeMesh: float = 0.05,
        minFacetAngle: float = 10,
        maxRadius: float = 0.05,
        path_data: str = ".",
        verbose: bool = False,
    ) -> bool:
        """
        Creates surfaces files of the TPMS and saves them in path_data directory
        """

        thickness = self.thickness * pi

        s_testplus = Generator(rve, thickness / 4.0, self.surface_function)
        s_testminus = Generator(rve, -thickness / 4.0, self.surface_function)
        s_plus = Generator(rve, thickness / 2.0, self.surface_function)
        s_minus = Generator(rve, -1.0 * thickness / 2.0, self.surface_function)

        print("mesh_surf_testplus  (1/4)", end="\r")
        mesh_surf_testplus = pygalmesh.generate_surface_mesh(
            s_testplus,
            min_facet_angle=minFacetAngle,
            max_radius_surface_delaunay_ball=maxRadius,
            max_facet_distance=sizeMesh,
            verbose=verbose,
        )

        print("mesh_surf_testminus (2/4)", end="\r")
        mesh_surf_testminus = pygalmesh.generate_surface_mesh(
            s_testminus,
            min_facet_angle=minFacetAngle,
            max_radius_surface_delaunay_ball=maxRadius,
            max_facet_distance=sizeMesh,
            verbose=verbose,
        )

        print("mesh_surf_plus      (4/4)", end="\r")
        mesh_surf_plus = pygalmesh.generate_surface_mesh(
            s_plus,
            min_facet_angle=minFacetAngle,
            max_radius_surface_delaunay_ball=maxRadius,
            max_facet_distance=sizeMesh,
            verbose=verbose,
        )

        print("mesh_surf_minus     (4/4)", end="\r")
        mesh_surf_minus = pygalmesh.generate_surface_mesh(
            s_minus,
            min_facet_angle=minFacetAngle,
            max_radius_surface_delaunay_ball=maxRadius,
            max_facet_distance=sizeMesh,
            verbose=verbose,
        )

        if not os.path.isdir(path_data):
            os.mkdir(path_data)
        mesh_surf_testplus.write(path_data + "/" + "tpms_testplus.stl")
        mesh_surf_testminus.write(path_data + "/" + "tpms_testminus.stl")
        mesh_surf_plus.write(path_data + "/" + "tpms_plus.stl")
        mesh_surf_minus.write(path_data + "/" + "tpms_minus.stl")

        return True

    def createTpms(self, path_data: str, rve: Rve) -> cq.Workplane:
        """
        Creates TPMS geometry (sheet or skeletal part) from surface
        files located in the directory given by path_data
        """

        if rve is None:
            raise ValueError("Please add an RVE to generate the TPMS")

        surf_tp = TopoDS_Shape()
        surf_tm = TopoDS_Shape()
        surf_p = TopoDS_Shape()
        surf_m = TopoDS_Shape()
        stl_reader = StlAPI_Reader()

        if not (
            stl_reader.Read(surf_tp, path_data + "/" + "tpms_testplus.stl")
            and stl_reader.Read(surf_tm, path_data + "/" + "tpms_testminus.stl")
            and stl_reader.Read(surf_p, path_data + "/" + "tpms_plus.stl")
            and stl_reader.Read(surf_m, path_data + "/" + "tpms_minus.stl")
        ):
            raise ValueError(
                "tpms_plus.stl, tpms_minus.stl, tpms_testplus.stl, tpms_testminus.stl files not found in '"
                + path_data
                + "' folder"
            )

        face_cut_tp = cq.Face(surf_tp)
        face_cut_tm = cq.Face(surf_tm)
        face_cut_p = cq.Face(surf_p)
        face_cut_m = cq.Face(surf_m)

        box = cq.Workplane("front").box(rve.dx, rve.dy, rve.dz)

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
        return cq.Workplane().add(return_object[0])


class Generator(pygalmesh.DomainBase):
    """
    Base class to generate TPMS geometry with a given surface function
    """
    def __init__(
        self,
        rve: Rve,
        height: float,
        surface_function: Callable[[float, float, float, float], float],
    ) -> None:
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = (
            math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2 + (0.5 * rve.dz) ** 2)
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
