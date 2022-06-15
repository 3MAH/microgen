"""
=============================================
TPMS (:mod:`microgen.shape.tpms`)
=============================================
"""
import math
import os
from typing import Callable, Union

import cadquery as cq
import pygalmesh
from numpy import cos, pi, sin
from OCP.StlAPI import StlAPI_Reader
from OCP.TopoDS import TopoDS_Shape

from ..operations import fuseParts, repeatGeometry, rescale
from ..rve import Rve
from .basicGeometry import BasicGeometry


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
        type_part: str = "sheet",
        thickness: float = 0,
        cell_size: Union[float, tuple[float, float, float], None] = None,
        repeat_cell: Union[int, tuple[int, int, int], None] = None,
        path_data: str = ".",
    ) -> None:
        """
        :param center: center of the geometry
        :param orientation: orientation of the geometry
        :param surface_function: tpms function or custom function (f(x, y, z, t) = 0)
        :param type_part: 'sheet' or 'skeletal'
        :param thickness: thickness of the tpms
        :param cell_size: By default, the tpms is generated for (1, 1, 1) dimensions but this can be modified by passing 'cell_size' scaling parameter (float or list of float for each dimension)
        :param repeat_cell: By default, the tpms is generated for one unit_cell. 'repeat_cell' parameter allows to repeat the geometry in the three dimensions
        :param path_data: folder where to store surfaces stl files required to generate the tpms
        """
        super().__init__(shape="TPMS", center=center, orientation=orientation)

        self.surface_function = surface_function

        if type_part != "sheet" and type_part != "skeletal":
            raise ValueError("type_part must be 'sheet' or 'skeletal'")
        self.type_part = type_part

        self.thickness = thickness

        if isinstance(cell_size, tuple) or cell_size is None:
            self.cell_size = cell_size
        elif cell_size is None:
            self.cell_size = (1, 1, 1)
        else:
            self.cell_size = (cell_size, cell_size, cell_size)

        if isinstance(repeat_cell, tuple) or repeat_cell is None:
            self.repeat_cell = repeat_cell
        else:
            self.repeat_cell = (repeat_cell, repeat_cell, repeat_cell)

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
        """
        Create TPMS surface for the corresponding isovalue and save file to path_data

        :param isovalue: height isovalue of the given tpms function
        :param filename: surface file name
        :param sizeMesh: max_facet_distance from `pygalmesh`_
        :param minFacetAngle: min_facet_angle from `pygalmesh`_
        :param maxRadius: max_radius_surface_delaunay_ball from `pygalmesh`_
        :param verbose: verbose from `pygalmesh`_

        .. _pygalmesh: https://github.com/meshpro/pygalmesh
        """
        surface = Generator(height=isovalue, surface_function=self.surface_function)
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
        Generates TPMS surface within a box for isovalue = 0

        :param sizeMesh: max_facet_distance from `pygalmesh`_
        :param minFacetAngle: min_facet_angle from `pygalmesh`_
        :param maxRadius: max_radius_surface_delaunay_ball from `pygalmesh`_
        :param verbose: verbose from `pygalmesh`_

        .. _pygalmesh: https://github.com/meshpro/pygalmesh
        """
        self.createSurface(
            isovalue=0,
            filename="surface.stl",
            sizeMesh=sizeMesh,
            minFacetAngle=minFacetAngle,
            maxRadius=maxRadius,
            verbose=verbose,
        )
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

        :param sizeMesh: max_facet_distance from `pygalmesh`_
        :param minFacetAngle: min_facet_angle from `pygalmesh`_
        :param maxRadius: max_radius_surface_delaunay_ball from `pygalmesh`_
        :param verbose: verbose from `pygalmesh`_

        .. _pygalmesh: https://github.com/meshpro/pygalmesh
        """

        thickness = self.thickness * pi
        self.createSurface(
            isovalue=thickness / 4.0,
            filename="tpms_testplus.stl",
            sizeMesh=sizeMesh,
            minFacetAngle=minFacetAngle,
            maxRadius=maxRadius,
            verbose=verbose,
        )
        self.createSurface(
            isovalue=-thickness / 4.0,
            filename="tpms_testminus.stl",
            sizeMesh=sizeMesh,
            minFacetAngle=minFacetAngle,
            maxRadius=maxRadius,
            verbose=verbose,
        )
        self.createSurface(
            isovalue=thickness / 2.0,
            filename="tpms_plus.stl",
            sizeMesh=sizeMesh,
            minFacetAngle=minFacetAngle,
            maxRadius=maxRadius,
            verbose=verbose,
        )
        self.createSurface(
            isovalue=-thickness / 2.0,
            filename="tpms_minus.stl",
            sizeMesh=sizeMesh,
            minFacetAngle=minFacetAngle,
            maxRadius=maxRadius,
            verbose=verbose,
        )

        surf_tp = TopoDS_Shape()
        surf_tm = TopoDS_Shape()
        surf_p = TopoDS_Shape()
        surf_m = TopoDS_Shape()
        stl_reader = StlAPI_Reader()
        stl_reader.Read(surf_tp, self.path_data + "/" + "tpms_testplus.stl")
        stl_reader.Read(surf_tm, self.path_data + "/" + "tpms_testminus.stl")
        stl_reader.Read(surf_p, self.path_data + "/" + "tpms_plus.stl")
        stl_reader.Read(surf_m, self.path_data + "/" + "tpms_minus.stl")

        face_cut_tp = cq.Face(surf_tp)
        face_cut_tm = cq.Face(surf_tm)
        face_cut_p = cq.Face(surf_p)
        face_cut_m = cq.Face(surf_m)

        box_wp = cq.Workplane("front").box(1, 1, 1)

        boxCut_wp = box_wp.split(face_cut_p)
        boxCut_wp = boxCut_wp.split(face_cut_m)

        boxWorkplanes = boxCut_wp.solids().all()  # type: list[cq.Workplane]

        listShapes = []  # type: list[tuple[int, cq.Shape]]

        for wp in boxWorkplanes:
            listShapes.append(
                (wp.split(face_cut_tp).split(face_cut_tm).solids().size(), wp.val())
            )

        sheet = [shape for (number, shape) in listShapes if number > 1]
        skeletal = [shape for (number, shape) in listShapes if number == 1]

        if self.type_part == "sheet":
            to_fuse = [cq.Shape(shape.wrapped) for shape in sheet]
            return_object = fuseParts(to_fuse, True)
        elif self.type_part == "skeletal":
            to_fuse = [cq.Shape(shape.wrapped) for shape in skeletal]
            return_object = fuseParts(to_fuse, False)

        if self.cell_size is not None:
            return_object = rescale(
                obj=return_object, scale=self.cell_size
            )
        if self.repeat_cell is not None:
            rve = Rve(*self.cell_size)
            return_object = repeatGeometry(
                unit_geom=return_object, rve=rve, grid=self.repeat_cell
            )

        return return_object.shape


class Generator(pygalmesh.DomainBase):
    """
    Base class to generate TPMS geometry with a given surface function
    Inherited from DomainBase class from `pygalmesh`_

    .. _pygalmesh: https://github.com/meshpro/pygalmesh
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
        self.waist_radius = math.sqrt(0.5**2 + 0.5**2)
        self.bounding_sphere_squared_radius = (
            math.sqrt(0.5**2 + 0.5**2 + 0.5**2) * 1.1
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
