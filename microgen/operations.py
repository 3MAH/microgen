"""
Boolean operations
"""

from typing import Union, Tuple, List, Sequence

import OCP
import cadquery as cq
import numpy as np
import pyvista as pv
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain

from .phase import Phase
from .rve import Rve


def _getRotationAxes(
    psi: float, theta: float, phi: float
) -> List[Tuple[float, float, float]]:
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
    obj: Union[cq.Shape, cq.Workplane],
    center: Union[np.ndarray, Tuple[float, float, float]],
    psi: float,
    theta: float,
    phi: float,
) -> Union[cq.Shape, cq.Workplane]:
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


def rescale(
    shape: cq.Shape, scale: Union[float, Tuple[float, float, float]]
) -> cq.Shape:
    """
    Rescale given object according to scale parameters [dim_x, dim_y, dim_z]

    :param shape: Shape
    :param scale: float or list of scale factor in each direction

    :return shape: rescaled Shape
    """
    return Phase.rescaleShape(shape, scale)


def fuseShapes(cqShapeList: List[cq.Shape], retain_edges: bool) -> cq.Shape:
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
            shape = upgrader.Shape()  # type: OCP.TopoDS_Shape            
            return cq.Shape(shape)            
        except Exception:
            return cq.Shape(occ_Solids)







def cutPhasesByShape(phaseList: List[Phase], cut_obj: cq.Shape) -> List[Phase]:
    """
    Cuts list of phases by a given shape

    :param phaseList: list of phases to cut
    :param cut_obj: cutting object

    :return phase_cut: final result
    """
    phase_cut: List[Phase] = []

    for phase in phaseList:
        cut = BRepAlgoAPI_Cut(phase.shape.wrapped, cut_obj.wrapped)
        if len(cq.Shape(cut.Shape()).Solids()) > 0:
            phase_cut.append(Phase(shape=cq.Shape(cut.Shape())))
    return phase_cut


def cutPhaseByShapeList(phaseToCut: Phase, cqShapeList: List[cq.Shape]) -> Phase:
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


def cutShapes(cqShapeList: List[cq.Shape], reverseOrder: bool = True) -> List[cq.Shape]:
    """
    Cuts list of shapes in the given order (or reverse) and fuse them.

    :param cqShapeList: list of CQ Shape to cut
    :param reverseOrder: bool, order for cutting shapes, when True: the last shape of the list is not cutted

    :return cutted_shapes: list of CQ Shape
    """
    cutted_shapes: List[cq.Shape] = []
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


def cutPhases(phaseList: List[Phase], reverseOrder: bool = True) -> List[Phase]:
    """
    Cuts list of shapes in the given order (or reverse) and fuse them.

    :param phaseList: list of phases to cut
    :param reverseOrder: bool, order for cutting shapes, when True: the last shape of the list is not cutted

    :return list of phases
    """
    shapeList = [phase.shape for phase in phaseList]
    cutted_shapes = cutShapes(shapeList, reverseOrder)

    return [Phase(shape=shape) for shape in cutted_shapes]


# def rasterShapeList(
#     cqShapeList: List[cq.Shape], rve: Rve, grid: List[int]
# ) -> Tuple[List[cq.Solid], List[List[cq.Solid]], List[float], List[cq.Vector]]:
#     """
#     Rasters shapes from shape list according to the rve divided by the given grid

#     :param cqShapeList: list of shapes to raster
#     :param rve: RVE divided by the given grid
#     :param grid: number of divisions in each direction [x, y, z]

#     :return flat_list: flatten list from occ_solids_list
#     :return occ_solids_list: list of list of solids
#     :return volume_list: list of shape volume (scalar)
#     :return center_list: list of shape center (vector)
#     """

#     occ_solids_list = []

#     for cqshape in cqShapeList:
#         wk_plane = cq.Workplane().add(cqshape.Solids())
#         xgrid = np.linspace(rve.x_min, rve.x_max, num=grid[0])
#         ygrid = np.linspace(rve.y_min, rve.y_max, num=grid[1])
#         zgrid = np.linspace(rve.z_min, rve.z_max, num=grid[2])
#         np.delete(xgrid, 0)
#         np.delete(ygrid, 0)
#         np.delete(zgrid, 0)
#         for i in xgrid:
#             Plane_x = cq.Face.makePlane(basePnt=(i, 0, 0), dir=(1, 0, 0))
#             wk_plane = wk_plane.split(cq.Workplane().add(Plane_x))
#         for j in ygrid:
#             Plane_y = cq.Face.makePlane(basePnt=(0, j, 0), dir=(0, 1, 0))
#             wk_plane = wk_plane.split(cq.Workplane().add(Plane_y))
#         for k in zgrid:
#             Plane_z = cq.Face.makePlane(basePnt=(0, 0, k), dir=(0, 0, 1))
#             wk_plane = wk_plane.split(cq.Workplane().add(Plane_z))

#         occ_solids_list.append(wk_plane.val().Solids())

#     flat_list = [item for sublist in occ_solids_list for item in sublist]
#     volume_list = [item.Volume() for sublist in occ_solids_list for item in sublist]
#     center_list = [item.Center() for sublist in occ_solids_list for item in sublist]
#     return (flat_list, occ_solids_list, volume_list, center_list)


def rasterPhase(
    phase: Phase, rve: Rve, grid: List[int], phasePerRaster: bool = True
) -> Union[Phase, List[Phase]]:
    """
    Rasters solids from phase according to the rve divided by the given grid

    :param phase: phase to raster
    :param rve: RVE divided by the given grid
    :param grid: number of divisions in each direction [x, y, z]
    :param phasePerRaster: if True, returns list of phases

    :return: Phase or list of Phases
    """
    solidList: List[cq.Solid] = phase.buildSolids(rve, grid)

    if phasePerRaster:
        return Phase.generatePhasePerRaster(solidList, rve, grid)
    return Phase(solids=solidList)


def repeatShape(unit_geom: cq.Shape, rve: Rve, grid: Tuple[int, int, int]) -> cq.Shape:
    """
    Repeats unit geometry in each direction according to the given grid

    :param unit_geom: Shape to repeat
    :param rve: RVE of the geometry to repeat
    :param grid: list of number of geometry repetitions in each direction

    :return: cq shape of the repeated geometry
    """
    return Phase.repeatShape(unit_geom, rve, grid)


def repeatPolyData(
    mesh: pv.PolyData, rve: Rve, grid: Tuple[int, int, int]
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
                        -rve.dim_x * (0.5 * grid[0] - 0.5 - i_x),
                        -rve.dim_y * (0.5 * grid[1] - 0.5 - i_y),
                        -rve.dim_z * (0.5 * grid[2] - 0.5 - i_z),
                    ),
                    inplace=True,
                )
                xyz_repeat.merge(new_mesh, inplace=True)
    return xyz_repeat
