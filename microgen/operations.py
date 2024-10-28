"""
Boolean operations
"""

from __future__ import annotations

from typing import Sequence

import cadquery as cq
import numpy as np
import OCP
import pyvista as pv
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain

from .phase import Phase
from .rve import Rve


def _getRotationAxes(
    psi: float, theta: float, phi: float
) -> list[tuple[float, float, float]]:
    """
    Retrieve the 3 Euler rotation axes

    :param psi: first Euler angle, in degrees
    :param theta: first Euler angle, in degrees
    :param phi: first Euler angle, in degrees

    :return: a list containing the three 3D Euler rotation axes
    """
    psi_rad, theta_rad, phi_rad = np.deg2rad((psi, theta, phi))
    return [
        (0.0, 0.0, 1.0),
        (np.cos(psi_rad), np.sin(psi_rad), 0.0),
        (
            np.sin(psi_rad) * np.sin(theta_rad),
            -np.sin(theta_rad) * np.cos(psi_rad),
            np.cos(theta_rad),
        ),
    ]


def rotateEuler(
    obj: cq.Shape | cq.Workplane,
    center: np.ndarray | tuple[float, float, float],
    psi: float,
    theta: float,
    phi: float,
) -> cq.Shape | cq.Workplane:
    """
    Rotates object according to XZX Euler angle convention

    :param obj: Object to rotate
    :param center: numpy array (x, y, z)
    :param psi: first Euler angle, in degrees
    :param theta: first Euler angle, in degrees
    :param phi: first Euler angle, in degrees

    :return: Rotated object
    """
    center_vector = cq.Vector(*center)
    z, u, z2 = _getRotationAxes(psi, theta, phi)
    for axis, angle in zip((z, u, z2), (psi, theta, phi)):
        obj = obj.rotate(center_vector, center_vector + cq.Vector(*axis), angle)
    return obj


def rotatePvEuler(
    obj: pv.PolyData,
    center: Sequence[float],
    psi: float,
    theta: float,
    phi: float,
) -> pv.PolyData:
    """
    Rotates object according to XZX Euler angle convention

    :param obj: Object to rotate
    :param center: numpy array (x, y, z)
    :param psi: first Euler angle, in degrees
    :param theta: second Euler angle, in degrees
    :param phi: third Euler angle, in degrees

    :return: Rotated object
    """
    z, u, z2 = _getRotationAxes(psi, theta, phi)

    rotated_obj = obj.rotate_vector(
        vector=z, angle=psi, point=tuple(center), inplace=False
    )
    rotated_obj.rotate_vector(vector=u, angle=theta, point=tuple(center), inplace=True)
    rotated_obj.rotate_vector(vector=z2, angle=phi, point=tuple(center), inplace=True)
    return rotated_obj


def rescale(shape: cq.Shape, scale: float | tuple[float, float, float]) -> cq.Shape:
    """
    Rescale given object according to scale parameters [dim_x, dim_y, dim_z]

    :param shape: Shape
    :param scale: float or list of scale factor in each direction

    :return shape: rescaled Shape
    """
    return Phase.rescaleShape(shape, scale)


def fuseShapes(cqShapeList: list[cq.Shape], retain_edges: bool) -> cq.Shape:
    """
    Fuse all shapes in cqShapeList

    :param cqShapeList: list of shapes to fuse
    :param retain_edges: retain intersecting edges

    :return fused object
    """

    occ_Solids = cqShapeList[0].wrapped
    for i in range(1, len(cqShapeList)):
        fuse = BRepAlgoAPI_Fuse(occ_Solids, cqShapeList[i].wrapped)
        occ_Solids = fuse.Shape()

    if retain_edges:
        return cq.Shape(occ_Solids)
    else:
        try:
            upgrader = ShapeUpgrade_UnifySameDomain(occ_Solids, True, True, True)
            upgrader.Build()
            shape: OCP.TopoDS_Shape = upgrader.Shape()
            return cq.Shape(shape)
        except Exception:
            return cq.Shape(occ_Solids)


def cutPhasesByShape(phaseList: list[Phase], cut_obj: cq.Shape) -> list[Phase]:
    """
    Cuts list of phases by a given shape

    :param phaseList: list of phases to cut
    :param cut_obj: cutting object

    :return phase_cut: final result
    """
    phase_cut: list[Phase] = []

    for phase in phaseList:
        cut = BRepAlgoAPI_Cut(phase.shape.wrapped, cut_obj.wrapped)
        if len(cq.Shape(cut.Shape()).Solids()) > 0:
            phase_cut.append(Phase(shape=cq.Shape(cut.Shape())))
    return phase_cut


def cutPhaseByShapeList(phaseToCut: Phase, cqShapeList: list[cq.Shape]) -> Phase:
    """
    Cuts a phase by a list of shapes

    :param phaseToCut: phase to cut
    :param cqShapeList: list of cutting shapes

    :return resultCut: cut phase
    """

    resultCut = phaseToCut.shape
    for shape in cqShapeList:
        cut = BRepAlgoAPI_Cut(resultCut.wrapped, shape.wrapped)
        resultCut = cq.Shape(cut.Shape())
    return Phase(shape=resultCut)


def cutShapes(cqShapeList: list[cq.Shape], reverseOrder: bool = True) -> list[cq.Shape]:
    """
    Cuts list of shapes in the given order (or reverse) and fuse them.

    :param cqShapeList: list of CQ Shape to cut
    :param reverseOrder: bool, order for cutting shapes, when True: the last shape of the list is not cut

    :return cutted_shapes: list of CQ Shape
    """
    cutted_shapes: list[cq.Shape] = []
    if reverseOrder:
        cqShapeList_inv = cqShapeList[::-1]
    else:
        cqShapeList_inv = cqShapeList

    cut_shape = cqShapeList_inv[0].copy()
    cutted_shapes.append(cut_shape)

    for shape in cqShapeList_inv[1::]:
        copy = shape.copy()
        cut = BRepAlgoAPI_Cut(copy.wrapped, cut_shape.wrapped)
        cutted_shapes.append(cq.Shape(cut.Shape()))

        fuse = BRepAlgoAPI_Fuse(cut_shape.wrapped, shape.wrapped)
        fused = fuse.Shape()
        upgrader = ShapeUpgrade_UnifySameDomain(fused, True, True, True)
        upgrader.Build()
        cut_shape = cq.Shape(upgrader.Shape())

    cutted_shapes.reverse()

    return cutted_shapes


def cutPhases(phaseList: list[Phase], reverseOrder: bool = True) -> list[Phase]:
    """
    Cuts list of shapes in the given order (or reverse) and fuse them.

    :param phaseList: list of phases to cut
    :param reverseOrder: bool, order for cutting shapes, when True: the last shape of the list is not cut

    :return list of phases
    """
    shapeList = [phase.shape for phase in phaseList]
    cutted_shapes = cutShapes(shapeList, reverseOrder)

    return [Phase(shape=shape) for shape in cutted_shapes]


def rasterPhase(
    phase: Phase, rve: Rve, grid: list[int], phasePerRaster: bool = True
) -> Phase | list[Phase]:
    """
    Rasters solids from phase according to the rve divided by the given grid

    :param phase: phase to raster
    :param rve: RVE divided by the given grid
    :param grid: number of divisions in each direction [x, y, z]
    :param phasePerRaster: if True, returns list of phases

    :return: Phase or list of Phases
    """
    solids: list[cq.Solid] = phase.split_solids(rve, grid)

    if phasePerRaster:
        return Phase.generate_phase_per_raster(solids, rve, grid)
    return Phase(solids=solids)


def repeat_shape(unit_geom: cq.Shape, rve: Rve, grid: tuple[int, int, int]) -> cq.Shape:
    """
    Repeats unit geometry in each direction according to the given grid

    :param unit_geom: Shape to repeat
    :param rve: RVE of the geometry to repeat
    :param grid: list of number of geometry repetitions in each direction

    :return: cq shape of the repeated geometry
    """
    return Phase.repeat_shape(unit_geom, rve, grid)


def repeatShape(unit_geom: cq.Shape, rve: Rve, grid: tuple[int, int, int]) -> cq.Shape:
    """
    Repeats unit geometry in each direction according to the given grid

    :param unit_geom: Shape to repeat
    :param rve: RVE of the geometry to repeat
    :param grid: list of number of geometry repetitions in each direction

    :return: cq shape of the repeated geometry
    """
    return Phase.repeatShape(unit_geom, rve, grid)


def repeatPolyData(
    mesh: pv.PolyData, rve: Rve, grid: tuple[int, int, int]
) -> pv.PolyData:
    """
    Repeats mesh in each direction according to the given grid

    :param mesh: pv.PolyData to repeat
    :param rve: RVE of the geometry to repeat
    :param grid: list of number of geometry repetitions in each direction

    :return: pv.PolyData of the repeated geometry
    """

    xyz_repeat = pv.PolyData()
    for i_x in range(grid[0]):
        for i_y in range(grid[1]):
            for i_z in range(grid[2]):
                new_mesh = mesh.copy()
                new_mesh.translate(
                    (
                        -rve.dim[0] * (0.5 * grid[0] - 0.5 - i_x),
                        -rve.dim[1] * (0.5 * grid[1] - 0.5 - i_y),
                        -rve.dim[2] * (0.5 * grid[2] - 0.5 - i_z),
                    ),
                    inplace=True,
                )
                xyz_repeat.merge(new_mesh, inplace=True)
    return xyz_repeat
