"""Boolean operations."""

from __future__ import annotations

import itertools
import warnings
from typing import TYPE_CHECKING, Sequence

import cadquery as cq
import numpy as np
import numpy.typing as npt
import pyvista as pv
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain

from .phase import Phase

if TYPE_CHECKING:
    import OCP
    from .rve import Rve


def _get_rotation_axes(
    psi: float,
    theta: float,
    phi: float,
) -> list[tuple[float, float, float]]:
    """Retrieve the 3 Euler rotation axes.

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


def rotate_euler(
    obj: cq.Shape | cq.Workplane,
    center: npt.NDArray[np.float64] | Sequence[float],
    angles: Sequence[float],
) -> cq.Shape | cq.Workplane:
    """Rotate object according to XZX Euler angle convention.

    :param obj: Object to rotate
    :param center: numpy array (x, y, z)
    :param angles: list of Euler angles (psi, theta, phi) in degrees

    :return: Rotated object
    """
    center_vector = cq.Vector(*center)
    psi, theta, phi = angles
    z, u, z2 = _get_rotation_axes(psi, theta, phi)
    for axis, angle in zip((z, u, z2), (psi, theta, phi)):
        obj = obj.rotate(center_vector, center_vector + cq.Vector(*axis), angle)
    return obj


def rotate_pv_euler(
    obj: pv.PolyData,
    center: Sequence[float],
    angles: Sequence[float],
) -> pv.PolyData:
    """Rotate object according to XZX Euler angle convention.

    :param obj: Object to rotate
    :param center: numpy array (x, y, z)
    :param angles: list of Euler angles (psi, theta, phi) in degrees

    :return: Rotated object
    """
    psi, theta, phi = angles
    z, u, z2 = _get_rotation_axes(psi, theta, phi)

    rotated_obj = obj.rotate_vector(
        vector=z,
        angle=psi,
        point=tuple(center),
        inplace=False,
    )
    rotated_obj.rotate_vector(vector=u, angle=theta, point=tuple(center), inplace=True)
    rotated_obj.rotate_vector(vector=z2, angle=phi, point=tuple(center), inplace=True)
    return rotated_obj


def rescale(shape: cq.Shape, scale: float | tuple[float, float, float]) -> cq.Shape:
    """Rescale given object according to scale parameters [dim_x, dim_y, dim_z].

    :param shape: Shape
    :param scale: float or list of scale factor in each direction

    :return shape: rescaled Shape
    """
    return Phase.rescaleShape(shape, scale)


def _unify_solids(solids: list[OCP.TopoDS_Shape]) -> cq.Shape:
    unify_edges = True
    unify_faces = True
    concat_bsplines = True
    upgrader = ShapeUpgrade_UnifySameDomain(
        solids,
        unify_edges,
        unify_faces,
        concat_bsplines,
    )
    upgrader.Build()
    shape: OCP.TopoDS_Shape = upgrader.Shape()
    return cq.Shape(shape)


def fuse_shapes(shapes: list[cq.Shape], *, retain_edges: bool) -> cq.Shape:
    """Fuse all shapes in cqShapeList.

    :param shapes: list of shapes to fuse
    :param retain_edges: retain intersecting edges

    :return fused object
    """
    fused = shapes[0].wrapped
    for i in range(1, len(shapes)):
        fused = BRepAlgoAPI_Fuse(fused, shapes[i].wrapped).Shape()

    if retain_edges:
        return cq.Shape(fused)

    try:
        return _unify_solids(fused)
    except Exception:  # which exception?
        return cq.Shape(fused)


def cut_phases_by_shape(phases: list[Phase], cut_obj: cq.Shape) -> list[Phase]:
    """Cut list of phases by a given shape.

    :param phases: list of phases to cut
    :param cut_obj: cutting object

    :return phase_cut: final result
    """
    phase_cut: list[Phase] = []

    for phase in phases:
        cut = cq.Shape(BRepAlgoAPI_Cut(phase.shape.wrapped, cut_obj.wrapped).Shape())
        if len(cut.Solids()) > 0:
            phase_cut.append(Phase(shape=cut))
    return phase_cut


def cut_phase_by_shape_list(phase_to_cut: Phase, shapes: list[cq.Shape]) -> Phase:
    """Cut a phase by a list of shapes.

    :param phase_to_cut: phase to cut
    :param shapes: list of cutting shapes

    :return resultCut: cut phase
    """
    result = phase_to_cut.shape
    for shape in shapes:
        result = cq.Shape(BRepAlgoAPI_Cut(result.wrapped, shape.wrapped).Shape())
    return Phase(shape=result)


def cut_shapes(shapes: list[cq.Shape], *, reverse_order: bool = True) -> list[cq.Shape]:
    """Cut list of shapes in the given order (or reverse) and fuse them.

    :param shapes: list of CQ Shape to cut
    :param reverse_order: bool, order for cutting shapes, \
        when True: the last shape of the list is not cut

    :return cutted_shapes: list of CQ Shape
    """
    cutted_shapes: list[cq.Shape] = []

    shapes_inv = reversed(shapes) if reverse_order else shapes

    cut_shape = shapes_inv[0].copy()
    cutted_shapes.append(cut_shape)

    for shape in shapes_inv[1::]:
        copy = shape.copy()
        cut = cq.Shape(BRepAlgoAPI_Cut(copy.wrapped, cut_shape.wrapped).Shape())
        cutted_shapes.append(cut)

        fused = BRepAlgoAPI_Fuse(cut_shape.wrapped, shape.wrapped).Shape()
        cut_shape = _unify_solids(fused)

    cutted_shapes.reverse()

    return cutted_shapes


def cut_phases(phases: list[Phase], *, reverse_order: bool = True) -> list[Phase]:
    """Cut list of shapes in the given order (or reverse) and fuse them.

    :param phases: list of phases to cut
    :param reverse_order: bool, order for cutting shapes, \
        when True: the last shape of the list is not cut

    :return list of phases
    """
    shapes = [phase.shape for phase in phases]
    cutted_shapes = cut_shapes(shapes, reverse_order=reverse_order)

    return [Phase(shape=shape) for shape in cutted_shapes]


def raster_phase(
    phase: Phase,
    rve: Rve,
    grid: list[int],
    *,
    phase_per_raster: bool = True,
) -> Phase | list[Phase]:
    """Raster solids from phase according to the rve divided by the given grid.

    :param phase: phase to raster
    :param rve: RVE divided by the given grid
    :param grid: number of divisions in each direction [x, y, z]
    :param phase_per_raster: if True, returns list of phases

    :return: Phase or list of Phases
    """
    solids: list[cq.Solid] = phase.split_solids(rve, grid)

    if phase_per_raster:
        return Phase.generate_phase_per_raster(solids, rve, grid)
    return Phase(solids=solids)


def repeat_shape(unit_geom: cq.Shape, rve: Rve, grid: tuple[int, int, int]) -> cq.Shape:
    """Repeat unit geometry in each direction according to the given grid.

    :param unit_geom: Shape to repeat
    :param rve: RVE of the geometry to repeat
    :param grid: list of number of geometry repetitions in each direction

    :return: cq shape of the repeated geometry
    """
    return Phase.repeat_shape(unit_geom, rve, grid)


def repeat_polydata(
    mesh: pv.PolyData,
    rve: Rve,
    grid: tuple[int, int, int],
) -> pv.PolyData:
    """Repeat mesh in each direction according to the given grid.

    :param mesh: pv.PolyData to repeat
    :param rve: RVE of the geometry to repeat
    :param grid: list of number of geometry repetitions in each direction

    :return: pv.PolyData of the repeated geometry
    """
    xyz_repeat = pv.PolyData()
    for idx in itertools.product(
        range(grid[0]),
        range(grid[1]),
        range(grid[2]),
    ):
        new_mesh = mesh.copy()
        xyz = [-rve.dim[i] * (0.5 * grid[i] - 0.5 - idx[i]) for i in range(3)]
        new_mesh.translate(xyz)
        xyz_repeat.merge(new_mesh)
    return xyz_repeat


# Deprecated functions
def rotateEuler(  # noqa: N802
    obj: cq.Shape | cq.Workplane,
    center: np.ndarray | tuple[float, float, float],
    psi: float,
    theta: float,
    phi: float,
) -> cq.Shape | cq.Workplane:
    """See rotate_euler.

    Deprecated in favor of rotate_euler.
    """
    warnings.warn(
        "rotateEuler is deprecated, use rotate_euler instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return rotate_euler(obj, center, (psi, theta, phi))


def rotatePvEuler(  # noqa: N802
    obj: pv.PolyData,
    center: Sequence[float],
    psi: float,
    theta: float,
    phi: float,
) -> pv.PolyData:
    """See rotatePvEuler.

    Deprecated in favor of rotatePvEuler.
    """
    warnings.warn(
        "rotatePvEuler is deprecated, use rotate_pv_euler instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return rotate_pv_euler(obj, center, (psi, theta, phi))


def fuseShapes(cqShapeList: list[cq.Shape], retain_edges: bool) -> cq.Shape:  # noqa: N802, N803, FBT001
    """See fuse_shapes.

    Deprecated in favor of fuse_shapes.
    """
    warnings.warn(
        "fuseShapes is deprecated, use fuse_shapes instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return fuse_shapes(cqShapeList, retain_edges=retain_edges)


def cutPhasesByShape(phaseList: list[Phase], cut_obj: cq.Shape) -> list[Phase]:  # noqa: N802, N803
    """See cut_phases_by_shape.

    Deprecated in favor of cut_phases_by_shape.
    """
    warnings.warn(
        "cutPhasesByShape is deprecated, use cut_phases_by_shape instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return cut_phases_by_shape(phaseList, cut_obj)


def cutPhaseByShapeList(phaseToCut: Phase, cqShapeList: list[cq.Shape]) -> Phase:  # noqa: N802, N803
    """See cut_phase_by_shape_list.

    Deprecated in favor of cut_phase_by_shape_list.
    """
    warnings.warn(
        "cutPhaseByShapeList is deprecated, use cut_phase_by_shape_list instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return cut_phase_by_shape_list(phaseToCut, cqShapeList)


def cutShapes(cqShapeList: list[cq.Shape], reverseOrder: bool = True) -> list[cq.Shape]:  # noqa: N802, N803, FBT001, FBT002
    """See cut_shapes.

    Deprecated in favor of cut_shapes.
    """
    warnings.warn(
        "cutShapes is deprecated, use cut_shapes instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return cut_shapes(cqShapeList, reverse_order=reverseOrder)


def cutPhases(phaseList: list[Phase], reverseOrder: bool = True) -> list[Phase]:  # noqa: N802, N803, FBT001, FBT002
    """See cut_phases.

    Deprecated in favor of cut_phases.
    """
    warnings.warn(
        "cutPhases is deprecated, use cut_phases instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return cut_phases(phaseList, reverse_order=reverseOrder)


def rasterPhase(  # noqa: N802
    phase: Phase,
    rve: Rve,
    grid: list[int],
    phasePerRaster: bool = True,  # noqa: N803, FBT001, FBT002
) -> Phase | list[Phase]:
    """See raster_phase.

    Deprecated in favor of raster_phase.
    """
    warnings.warn(
        "rasterPhase is deprecated, use raster_phase instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return raster_phase(phase, rve, grid, phase_per_raster=phasePerRaster)


def repeatShape(unit_geom: cq.Shape, rve: Rve, grid: tuple[int, int, int]) -> cq.Shape:  # noqa: N802
    """See repeat_shape.

    Deprecated in favor of repeat_shape.
    """
    warnings.warn(
        "repeatShape is deprecated, use repeat_shape instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return repeat_shape(unit_geom, rve, grid)


def repeatPolyData(  # noqa: N802
    mesh: pv.PolyData,
    rve: Rve,
    grid: tuple[int, int, int],
) -> pv.PolyData:
    """See repeat_polydata.

    Deprecated in favor of repeat_polydata.
    """
    warnings.warn(
        "repeatPolyData is deprecated, use repeat_polydata instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return repeat_polydata(mesh, rve, grid)
