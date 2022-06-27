"""
=============================================
TPMS (:mod:`microgen.shape.tpms`)
=============================================
"""
import math
import os
from typing import Callable, Union

import numpy as np
import cadquery as cq
import pyvista as pv
from numpy import cos, pi, sin
from OCP.StlAPI import StlAPI_Reader
from OCP.TopoDS import TopoDS_Shape

from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_MakeVertex,
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakePolygon,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_Sewing,
    BRepBuilderAPI_Copy,
    BRepBuilderAPI_GTransform,
    BRepBuilderAPI_Transform,
    BRepBuilderAPI_Transformed,
    BRepBuilderAPI_RightCorner,
    BRepBuilderAPI_RoundCorner,
    BRepBuilderAPI_MakeSolid,
)

from ..operations import fuseParts, rescale, repeatGeometry, repeatShape
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
        cell_size: Union[float, tuple[float, float, float], None] = 1,
        repeat_cell: Union[int, tuple[int, int, int], None] = None
    ) -> None:
        """
        :param center: center of the geometry
        :param orientation: orientation of the geometry
        :param surface_function: tpms function or custom function (f(x, y, z, t) = 0)
        :param type_part: 'sheet' or 'skeletal'
        :param thickness: thickness of the tpms
        :param cell_size: By default, the tpms is generated for (1, 1, 1) dimensions but this can be modified by passing 'cell_size' scaling parameter (float or list of float for each dimension)
        :param repeat_cell: By default, the tpms is generated for one unit_cell. 'repeat_cell' parameter allows to repeat the geometry in the three dimensions
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
        
    def createSurface(
        self,
        isovalue: float = 0,
        nSample: int = 20,
        smoothing: int = 100,
        verbose: bool = False,
    ) -> cq.Shell:
        """
        Create TPMS surface for the corresponding isovalue

        :param isovalue: height isovalue of the given tpms function
        :param nSample: surface file name
        :param smoothing: smoothing loop iterations
        """
        x_min, y_min, z_min = -0.5, -0.5, -0.5
        grid = pv.UniformGrid(
            dims=(nSample, nSample, nSample),
            spacing=(1./(nSample-1) , 1./(nSample-1), 1./(nSample-1)),
            origin=(x_min, y_min, z_min),
        )
        x, y, z = grid.points.T
        
        surface_function=self.surface_function
        surface = surface_function(x,y,z)
        mesh = grid.contour(1, scalars=surface, method='flying_edges', rng=(isovalue,1))
        if smoothing > 0:
            mesh = mesh.smooth(n_iter=smoothing)
        list_of_Triangles = mesh.faces.reshape(-1, 4)[:,1:]
        list_of_Triangles = np.c_[list_of_Triangles,list_of_Triangles[:,0]]

        faces = []
        for ixs in list_of_Triangles:
            lines = []
            for v1, v2 in zip(ixs, ixs[1:]):
                vertice_coords1 = mesh.points[v1]
                vertice_coords2 = mesh.points[v2]
                lines.append(
                    cq.Edge.makeLine(
                        cq.Vector(*vertice_coords1), cq.Vector(*vertice_coords2)
                    )
                )
            wire = cq.Wire.assembleEdges(lines)
            faces.append(cq.Face.makeFromWires(wire))

        return cq.Shell.makeShell(faces)

    def createSurfaces(
        self,
        isovalues: list[float] = [0],
        nSample: int = 20,
        smoothing: int = 100,
        verbose: bool = False,
    ) -> list[cq.Shell]:
        """
        Create TPMS surface for the corresponding isovalue
        
        :param numsber_surfaces: number of surfaces
        :param isovalues: height isovalues of the given tpms function
        :param nSample: surface file name
        :param smoothing: smoothing loop iterations
        """
        x_min, y_min, z_min = -0.5, -0.5, -0.5
        grid = pv.UniformGrid(
            dims=(nSample, nSample, nSample),
            spacing=(1./(nSample-1) , 1./(nSample-1), 1./(nSample-1)),
            origin=(x_min, y_min, z_min),
        )
        x, y, z = grid.points.T
        
        surface_function=self.surface_function
        surface = surface_function(x,y,z)

        shells = []
        for isovalue in isovalues:
            mesh = grid.contour(1, scalars=surface, method='flying_edges', rng=(isovalue,1))
            if smoothing > 0:
                mesh = mesh.smooth(n_iter=smoothing)
            list_of_Triangles = mesh.faces.reshape(-1, 4)[:,1:]
            list_of_Triangles = np.c_[list_of_Triangles,list_of_Triangles[:,0]]

            faces = []
            for ixs in list_of_Triangles:
                lines = []
                for v1, v2 in zip(ixs, ixs[1:]):
                    vertice_coords1 = mesh.points[v1]
                    vertice_coords2 = mesh.points[v2]
                    lines.append(
                        cq.Edge.makeLine(
                            cq.Vector(*vertice_coords1), cq.Vector(*vertice_coords2)
                        )
                    )
                wire = cq.Wire.assembleEdges(lines)
                faces.append(cq.Face.makeFromWires(wire))
            shells.append(cq.Shell.makeShell(faces))
            
        return shells

    def generateSurface(
        self,
        isovalue: float = 0.,
        nSample: int = 20,
        smoothing: int = 100,
        verbose: bool = False,
    ) -> cq.Shape:

        shell = self.createSurface(isovalue=isovalue, nSample=nSample, smoothing=smoothing)
        return_object = cq.Shape(shell.wrapped)
        if self.cell_size != (1.0, 1.0, 1.0):
            transform_mat = cq.Matrix(
                [
                    [self.cell_size[0], 0, 0, 0],
                    [0, self.cell_size[1], 0, 0],
                    [0, 0, self.cell_size[2], 0],
                ]
            )
            return_object = return_object.transformGeometry(transform_mat)
        
        if self.repeat_cell != (1, 1, 1):
            return_object = repeatShape(
                unit_geom=return_object, rve=Rve(*self.cell_size), grid=self.repeat_cell
            )
        return return_object

    def generate(
        self,
        nSample: int = 20,
        smoothing: int = 100,
        verbose: bool = False,
    ) -> cq.Shape:
        """
        Creates thick TPMS geometry (sheet or skeletal part) from surface

        :param isovalue: height isovalue of the given tpms function
        :param nSample: surface file name
        :param smoothing: smoothing loop iterations
        """

        isovalues= [-self.thickness,-self.thickness/3., self.thickness/3., self.thickness]
        shells = self.createSurfaces(isovalues=isovalues, nSample=nSample, smoothing=smoothing)

        face_cut_tp = shells[2]
        face_cut_tm = shells[1]
        face_cut_p = shells[3]
        face_cut_m = shells[0]

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


# Lidinoid -> 0.5*(sin(2*x)*cos(y)*sin(z) + sin(2*y)*cos(z)*sin(x) + sin(2*z)*cos(x)*sin(y)) - 0.5*(cos(2*x)*cos(2*y) + cos(2*y)*cos(2*z) + cos(2*z)*cos(2*x)) + 0.15 = 0


def gyroid(x: float, y: float, z: float) -> float:
    """
    :math:`sin(2 \pi x) cos(2 \pi y) + sin(2 \pi y) cos(2 \pi z) + sin(2 \pi z) cos(2 \pi x) = 0`
    """
    return (
        sin(2 * pi * x) * cos(2 * pi * y)
            + sin(2 * pi * y) * cos(2 * pi * z)
            + sin(2 * pi * z) * cos(2 * pi * x)
    )

def schwarzP(x: float, y: float, z: float, height: float) -> float:
    """
    :math:`cos(2 \pi x) + cos(2 \pi y) + cos(2 \pi z) = 0`
    """
    return cos(2 * pi * x) + cos(2 * pi * y) + cos(2 * pi * z)


def schwarzD(x: float, y: float, z: float) -> float:
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
    return a + b + c + d


def neovius(x: float, y: float, z: float) -> float:
    """
    :math:`3 cos(2 \pi x) + cos(2 \pi y) + cos(2 \pi z) + \
           4 cos(2 \pi x) cos(2 \pi y) cos(2 \pi z) = 0`
    """
    a = 3 * cos(2 * pi * x) + cos(2 * pi * y) + cos(2 * pi * z)
    b = 4 * cos(2 * pi * x) * cos(2 * pi * y) * cos(2 * pi * z)

    return a + b


def schoenIWP(x: float, y: float, z: float) -> float:
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

    return a - b


def schoenFRD(x: float, y: float, z: float) -> float:
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
    return a - b


def fischerKochS(x: float, y: float, z: float) -> float:
    """
    :math:`cos(4 \pi x) sin(2 \pi y) cos(2 \pi z) + \
           cos(2 \pi x) cos(4 \pi y) sin(2 \pi z) + \
           sin(2 \pi x) cos(2 \pi y) cos(4 \pi z) = 0`
    """
    a = cos(4 * pi * x) * sin(2 * pi * y) * cos(2 * pi * z)
    b = cos(2 * pi * x) * cos(4 * pi * y) * sin(2 * pi * z)
    c = sin(2 * pi * x) * cos(2 * pi * y) * cos(4 * pi * z)

    return a + b + c


def pmy(x: float, y: float, z: float) -> float:
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

    return a + b + c + d


def honeycomb(x: float, y: float, z: float) -> float:
    """
    :math:`sin(2 \pi x) cos(2 \pi y) + sin(2 \pi y) + cos(2 \pi z) = 0`
    """
    return (
        sin(2 * pi * x) * cos(2 * pi * y) + sin(2 * pi * y) + cos(2 * pi * z)
    )
