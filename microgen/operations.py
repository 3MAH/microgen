"""
Boolean operations
"""

import os
from typing import Union

import cadquery as cq
import numpy as np
import OCP
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain

from .phase import Phase
from .rve import Rve


def rotateEuler(
    object: Union[cq.Shape, cq.Workplane],
    center: Union[np.ndarray, tuple[float, float, float]],
    psi: float,
    theta: float,
    phi: float,
) -> Union[cq.Shape, cq.Workplane]:
    """
    Rotates object according to XZX Euler angle convention

    :param object: Object to rotate
    :param center: numpy array (x, y, z)
    :param psi, theta, phi: Euler angles

    :return object_r: Rotated object
    """

    u = np.array([0.0, 0.0, 1.0])
    u = np.array([np.cos(psi * np.pi / 180.0), np.sin(psi * np.pi / 180.0), 0.0])
    z2 = np.array(
        [
            np.sin(psi * np.pi / 180.0) * np.sin(theta * np.pi / 180.0),
            -np.sin(theta * np.pi / 180.0) * np.cos(psi * np.pi / 180.0),
            np.cos(theta * np.pi / 180.0),
        ]
    )

    object_r = object.rotate(
        cq.Vector(center[0], center[1], center[2]),
        cq.Vector(center[0], center[1], center[2] + 1.0),
        psi,
    )
    object_r = object_r.rotate(
        cq.Vector(center[0], center[1], center[2]),
        cq.Vector(center[0] + u[0], center[1] + u[1], center[2] + u[2]),
        theta,
    )
    object_r = object_r.rotate(
        cq.Vector(center[0], center[1], center[2]),
        cq.Vector(center[0] + z2[0], center[1] + z2[1], center[2] + z2[2]),
        phi,
    )
    return object_r


def rescale(
    obj: Phase, scale: Union[float, tuple[float, float, float]]
) -> Phase:
    """
    Rescale given object according to scale parameters [dim_x, dim_y, dim_z]

    :param obj: Phase
    :param scale: float or list of scale factor in each direction

    :return shape: rescaled Phase
    """
    if isinstance(scale, float):
        scale = (scale, scale, scale)

    center = obj.centerOfMass

    # move the shape at (0, 0, 0) to rescale it
    shape = obj.shape.move(cq.Location(cq.Vector(-center[0], -center[1], -center[2])))

    # then move it back to its center with transform Matrix
    transform_mat = cq.Matrix(
        [
            [scale[0], 0, 0, center[0]],
            [0, scale[1], 0, center[1]],
            [0, 0, scale[2], center[2]],
        ]
    )
    shape = shape.transformGeometry(transform_mat)

    return Phase(shape=shape)


def removeEmptyLines(filename: str) -> None:
    """
    Removes empty lines of the given file

    :param filename: file where to remove empty lines
    """
    if not os.path.isfile(filename):
        print("{} does not exist ".format(filename))
        return
    with open(filename) as filehandle:
        lines = filehandle.readlines()

    with open(filename, "w") as filehandle:
        lines = filter(lambda x: x.strip(), lines)
        filehandle.writelines(lines)


def fuseParts(cqShapeList: list[cq.Shape], retain_edges: bool) -> Phase:
    """
    Fuse all shapes in cqShapeList

    :param cqShapeList: list of shapes to fuse
    :param retain_edges: retain intersecting edges

    Returns
    -------
    cq.Shape : fused object
    occ_solids_list : list of list of solids
    """

    occ_solids_list = [s.Solids() for s in cqShapeList]

    occ_Solids = cqShapeList[0].wrapped
    for i in range(1, len(cqShapeList)):
        fuse = BRepAlgoAPI_Fuse(occ_Solids, cqShapeList[i].wrapped)
        occ_Solids = fuse.Shape()

    if retain_edges:
        phase = Phase(shape=cq.Shape(occ_Solids), solids=occ_solids_list)
        return phase
    else:
        upgrader = ShapeUpgrade_UnifySameDomain(occ_Solids, True, True, True)
        upgrader.Build()
        shape = upgrader.Shape()  # type: OCP.TopoDS_Shape
        return Phase(shape=cq.Shape(shape), solids=[[cq.Solid(shape)]])


def cutPhasesByShape(phaseList: list[Phase], cut_obj: cq.Shape) -> list[Phase]:
    """
    Cuts list of phases by a given shape

    :param phaseList: list of phases to cut
    :param cut_obj: cutting object

    :return phase_cut: final result
    """
    phase_cut = []  # type: list[Phase]

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

    :return resultCut: cutted phase
    """

    resultCut = phaseToCut.shape
    for shape in cqShapeList:
        cut = BRepAlgoAPI_Cut(resultCut.wrapped, shape.wrapped)
        resultCut = cq.Shape(cut.Shape())
    return Phase(shape=resultCut)


def cutParts(cqShapeList: list[cq.Shape], reverseOrder: bool = True) -> list[Phase]:
    """
    Cuts list of shapes in the given order (or reverse) and fuse them.

    :param cqShapeList: list of CQ Shape to cut
    :param reverseOrder: bool, order for cutting shapes, when True: the last shape of the list is not cutted

    :return phase_cut: list of phases
    """
    phase_cut = []  # type: list[Phase]
    if reverseOrder:
        cqShapeList_inv = cqShapeList[::-1]
    else:
        cqShapeList_inv = cqShapeList

    cut_obj = cqShapeList_inv[0].copy()
    phase_cut.append(Phase(shape=cut_obj))

    for shape in cqShapeList_inv[1::]:
        copy = shape.copy()
        cut = BRepAlgoAPI_Cut(copy.wrapped, cut_obj.wrapped)
        phase_cut.append(Phase(shape=cq.Shape(cut.Shape())))

        fuse = BRepAlgoAPI_Fuse(cut_obj.wrapped, shape.wrapped)
        fused = fuse.Shape()
        upgrader = ShapeUpgrade_UnifySameDomain(fused, True, True, True)
        upgrader.Build()
        cut_obj = cq.Shape(upgrader.Shape())

    phase_cut.reverse()

    return phase_cut


def rasterShapeList(
    cqShapeList: list[cq.Shape], rve: Rve, grid: list[int]
) -> tuple[list[cq.Solid], list[list[cq.Solid]], list[float], list[cq.Vector]]:
    """
    Rasters shapes from shape list according to the rve divided by the given grid

    :param cqShapeList: list of shapes to raster
    :param rve: RVE divided by the given grid
    :param grid: number of divisions in each direction [x, y, z]

    :return flat_list: flatten list from occ_solids_list
    :return occ_solids_list: list of list of solids
    :return volume_list: list of shape volume (scalar)
    :return center_list: list of shape center (vector)
    """

    occ_solids_list = []

    for cqshape in cqShapeList:
        wk_plane = cq.Workplane().add(cqshape.Solids())
        xgrid = np.linspace(0.0, rve.dx, num=grid[0])
        ygrid = np.linspace(0.0, rve.dy, num=grid[1])
        zgrid = np.linspace(0.0, rve.dz, num=grid[2])
        np.delete(xgrid, 0)
        np.delete(ygrid, 0)
        np.delete(zgrid, 0)
        for i in xgrid:
            Plane_x = cq.Face.makePlane(basePnt=(i, 0, 0), dir=(1, 0, 0))
            wk_plane = wk_plane.split(cq.Workplane().add(Plane_x))
        for j in ygrid:
            Plane_y = cq.Face.makePlane(basePnt=(0, j, 0), dir=(0, 1, 0))
            wk_plane = wk_plane.split(cq.Workplane().add(Plane_y))
        for k in zgrid:
            Plane_z = cq.Face.makePlane(basePnt=(0, 0, k), dir=(0, 0, 1))
            wk_plane = wk_plane.split(cq.Workplane().add(Plane_z))

        occ_solids_list.append(wk_plane.val().Solids())

    flat_list = [item for sublist in occ_solids_list for item in sublist]
    volume_list = [item.Volume() for sublist in occ_solids_list for item in sublist]
    center_list = [item.Center() for sublist in occ_solids_list for item in sublist]
    return (flat_list, occ_solids_list, volume_list, center_list)


def rasterPhase(
    phase: Phase, rve: Rve, grid: list[int], phasePerRaster: bool = True
) -> Union[Phase, list[Phase]]:
    """
    Rasters solids from phase according to the rve divided by the given grid

    :param phase: phase to raster
    :param rve: RVE divided by the given grid
    :param grid: number of divisions in each direction [x, y, z]

    :return: Phase
    """
    solidList = []  # type: list[cq.Solid]

    for solid in phase.solids:
        wk_plane = cq.Workplane().add(solid)
        xgrid = np.linspace(0.0, rve.dx, num=grid[0])
        ygrid = np.linspace(0.0, rve.dy, num=grid[1])
        zgrid = np.linspace(0.0, rve.dz, num=grid[2])
        np.delete(xgrid, 0)
        np.delete(ygrid, 0)
        np.delete(zgrid, 0)
        for i in xgrid:
            Plane_x = cq.Face.makePlane(basePnt=(i, 0, 0), dir=(1, 0, 0))
            wk_plane = wk_plane.split(cq.Workplane().add(Plane_x))
        for j in ygrid:
            Plane_y = cq.Face.makePlane(basePnt=(0, j, 0), dir=(0, 1, 0))
            wk_plane = wk_plane.split(cq.Workplane().add(Plane_y))
        for k in zgrid:
            Plane_z = cq.Face.makePlane(basePnt=(0, 0, k), dir=(0, 0, 1))
            wk_plane = wk_plane.split(cq.Workplane().add(Plane_z))

        for subsolid in wk_plane.val().Solids():
            solidList.append(subsolid)

    if phasePerRaster:
        solids_phases = [[] for _ in range(grid[0] * grid[1] * grid[2])]  # type: list[list[cq.Solid]]
        for solid in solidList:
            center = solid.Center()
            i = int(round(center.x / (rve.dx / grid[0])))
            j = int(round(center.y / (rve.dy / grid[1])))
            k = int(round(center.z / (rve.dz / grid[2])))
            ind = i + grid[0] * j + grid[0] * grid[1] * k
            solids_phases[ind].append(solid)
        return [Phase(solids=solids) for solids in solids_phases if len(solids) > 0]
    else:
        return Phase(solids=solidList)


def repeatGeometry(unit_geom: Phase, rve: Rve, grid: tuple[int, int, int]) -> Phase:
    """
    Repeats unit geometry in each direction according to the given grid

    :param unit_geom: Geometry to repeat
    :param rve: RVE of the geometry to repeat
    :param grid: list of number of geometry repetitions in each direction

    :return CQ_Compound: cq compound of the repeated geometry
    """

    xyz_repeat = cq.Assembly()
    for i_x in range(grid[0]):
        for i_y in range(grid[1]):
            for i_z in range(grid[2]):
                xyz_repeat.add(
                    unit_geom.shape,
                    loc=cq.Location(
                        cq.Vector(i_x * rve.dim_x, i_y * rve.dim_y, i_z * rve.dim_z)
                    ),
                )
    compound = xyz_repeat.toCompound()
    return Phase(shape=cq.Shape(compound.wrapped))
