"""Boolean operations.

All OCP imports are deferred so the module loads cleanly in a mesh-only
install.  Functions that use OCCT raise a clear ``ImportError`` (via
:func:`microgen.cad.require_cad`) if the ``[cad]`` extra isn't installed.
"""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING, overload

import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.spatial.transform import Rotation

from .cad import CadShape, require_cad

if TYPE_CHECKING:
    from collections.abc import Sequence

    from OCP.TopoDS import TopoDS_Shape

    from .phase import Phase
    from .rve import Rve

    Rotatable = CadShape | pv.PolyData


@overload
def rotate(
    obj: CadShape,
    center: npt.NDArray[np.float64] | Sequence[float],
    rotation: Rotation,
) -> CadShape: ...


@overload
def rotate(
    obj: pv.PolyData,
    center: npt.NDArray[np.float64] | Sequence[float],
    rotation: Rotation,
) -> pv.PolyData: ...


def rotate(
    obj: Rotatable,
    center: npt.NDArray[np.float64] | Sequence[float],
    rotation: Rotation,
) -> Rotatable:
    """Rotate object according to given rotation.

    Supports :class:`microgen.cad.CadShape` and :class:`pyvista.PolyData`.

    :param obj: Object to rotate
    :param center: numpy array (x, y, z)
    :param rotation: scipy Rotation object

    :return: Rotated object (same type as input)

    :raises ValueError: if object type is not supported

    """
    rotvec = rotation.as_rotvec(degrees=True)
    angle = np.linalg.norm(rotvec)
    if angle == 0:
        return obj
    axis = rotvec / angle

    if isinstance(obj, CadShape):
        return obj.rotate(center, axis, float(angle))
    if isinstance(obj, pv.PolyData):
        return obj.rotate_vector(axis, angle, point=tuple(center))

    err_msg = f"rotate(): object type {type(obj).__name__} not supported."
    raise ValueError(err_msg)


@overload
def rotate_euler(
    obj: CadShape,
    center: npt.NDArray[np.float64] | Sequence[float],
    angles_or_rotation: Sequence[float] | Rotation,
) -> CadShape: ...


@overload
def rotate_euler(
    obj: pv.PolyData,
    center: npt.NDArray[np.float64] | Sequence[float],
    angles_or_rotation: Sequence[float] | Rotation,
) -> pv.PolyData: ...


def rotate_euler(
    obj: Rotatable,
    center: npt.NDArray[np.float64] | Sequence[float],
    angles_or_rotation: Sequence[float] | Rotation,
) -> Rotatable:
    """Rotate object according to ZXZ Euler angle convention.

    Accepts :class:`~microgen.cad.CadShape` or :class:`pyvista.PolyData`.

    :param obj: Object to rotate
    :param center: numpy array (x, y, z)
    :param angles_or_rotation: list of Euler angles (psi, theta, phi) in
        degrees, or a scipy ``Rotation`` object

    :return: Rotated object
    """
    if isinstance(angles_or_rotation, Rotation):
        rotation = angles_or_rotation
    else:
        rotation = Rotation.from_euler("ZXZ", angles_or_rotation, degrees=True)
    return rotate(obj, center, rotation)


def rotate_pv_euler(
    obj: pv.PolyData,
    center: Sequence[float],
    angles_or_rotation: Sequence[float] | Rotation,
) -> pv.PolyData:
    """Rotate object according to ZXZ Euler angle convention.

    :param obj: Object to rotate
    :param center: numpy array (x, y, z)
    :param angles_or_rotation: list of Euler angles (psi, theta, phi) in
        degrees, or a scipy ``Rotation`` object

    :return: Rotated object
    """
    if isinstance(angles_or_rotation, Rotation):
        rotation = angles_or_rotation
    else:
        rotation = Rotation.from_euler("ZXZ", angles_or_rotation, degrees=True)

    rotvec = rotation.as_rotvec(degrees=True)
    angle = np.linalg.norm(rotvec)
    if angle == 0:
        return obj
    axis = rotvec / angle

    return obj.rotate_vector(
        vector=axis,
        angle=angle,
        point=tuple(center),
        inplace=False,
    )


def rescale(shape: CadShape, scale: float | tuple[float, float, float]) -> CadShape:
    """Rescale given object according to scale parameters [dim_x, dim_y, dim_z]."""
    from .phase import Phase  # noqa: PLC0415

    return Phase.rescale_shape(shape, scale)


def _unify_solids(shape: TopoDS_Shape) -> CadShape:
    require_cad()
    from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain  # noqa: PLC0415

    unify_edges = True
    unify_faces = True
    concat_bsplines = True
    upgrader = ShapeUpgrade_UnifySameDomain(
        shape,
        unify_edges,
        unify_faces,
        concat_bsplines,
    )
    upgrader.Build()
    return CadShape(upgrader.Shape())


def fuse_shapes(shapes: list[CadShape], *, retain_edges: bool) -> CadShape:
    """Fuse all shapes in the list.

    :param shapes: list of shapes to fuse
    :param retain_edges: retain intersecting edges

    :return: fused shape
    """
    require_cad()
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse  # noqa: PLC0415

    fused = shapes[0].wrapped
    for i in range(1, len(shapes)):
        fused = BRepAlgoAPI_Fuse(fused, shapes[i].wrapped).Shape()

    if retain_edges:
        return CadShape(fused)

    try:
        return _unify_solids(fused)
    except Exception:  # noqa: BLE001
        return CadShape(fused)


def cut_phases_by_shape(phases: list[Phase], cut_obj: CadShape) -> list[Phase]:
    """Cut list of phases by a given shape.

    :param phases: list of phases to cut
    :param cut_obj: cutting object

    :return phase_cut: final result
    """
    require_cad()
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut  # noqa: PLC0415

    from .phase import Phase  # noqa: PLC0415

    phase_cut: list[Phase] = []

    for phase in phases:
        if phase.shape is None:
            continue
        cut = CadShape(BRepAlgoAPI_Cut(phase.shape.wrapped, cut_obj.wrapped).Shape())
        if len(cut.solids()) > 0:
            phase_cut.append(Phase(shape=cut))
    return phase_cut


def cut_phase_by_shape_list(phase_to_cut: Phase, shapes: list[CadShape]) -> Phase:
    """Cut a phase by a list of shapes.

    :param phase_to_cut: phase to cut
    :param shapes: list of cutting shapes

    :return resultCut: cut phase
    """
    require_cad()
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut  # noqa: PLC0415

    from .phase import Phase  # noqa: PLC0415

    result = phase_to_cut.shape
    if result is None:
        err_msg = "phase_to_cut has no shape to cut"
        raise ValueError(err_msg)
    for shape in shapes:
        result = CadShape(BRepAlgoAPI_Cut(result.wrapped, shape.wrapped).Shape())
    return Phase(shape=result)


def cut_shapes(shapes: list[CadShape], *, reverse_order: bool = True) -> list[CadShape]:
    """Cut list of shapes in the given order (or reverse) and fuse them.

    :param shapes: list of shapes to cut
    :param reverse_order: bool, order for cutting shapes, \
        when True: the last shape of the list is not cut

    :return cutted_shapes: list of CadShape
    """
    require_cad()
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse  # noqa: PLC0415

    cutted_shapes: list[CadShape] = []

    shapes_inv = list(reversed(shapes)) if reverse_order else list(shapes)

    cut_shape = shapes_inv[0].copy()
    cutted_shapes.append(cut_shape)

    for shape in shapes_inv[1:]:
        copy_shape = shape.copy()
        cut = CadShape(BRepAlgoAPI_Cut(copy_shape.wrapped, cut_shape.wrapped).Shape())
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
    from .phase import Phase  # noqa: PLC0415

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
    from .phase import Phase  # noqa: PLC0415

    solids = phase.split_solids(rve, grid)

    if phase_per_raster:
        return Phase.generate_phase_per_raster(solids, rve, grid)
    return Phase(solids=solids)


def repeat_shape(unit_geom: CadShape, rve: Rve, grid: tuple[int, int, int]) -> CadShape:
    """Repeat unit geometry in each direction according to the given grid.

    :param unit_geom: Shape to repeat
    :param rve: RVE of the geometry to repeat
    :param grid: list of number of geometry repetitions in each direction

    :return: cq shape of the repeated geometry
    """
    from .phase import Phase  # noqa: PLC0415

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
