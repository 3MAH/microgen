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

    Rotatable = CadShape | pv.DataSet


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


@overload
def rotate(
    obj: pv.UnstructuredGrid,
    center: npt.NDArray[np.float64] | Sequence[float],
    rotation: Rotation,
) -> pv.UnstructuredGrid: ...


def rotate(
    obj: Rotatable,
    center: npt.NDArray[np.float64] | Sequence[float],
    rotation: Rotation,
) -> Rotatable:
    """Rotate object according to given rotation.

    Supports :class:`microgen.cad.CadShape` and any :class:`pyvista.DataSet`
    (``PolyData``, ``UnstructuredGrid``, ``StructuredGrid``, …) — i.e. any
    pyvista mesh exposing ``rotate_vector``.

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
    if isinstance(obj, pv.DataSet):
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


def rescale(shape: CadShape, scale: float | tuple[float, float, float]) -> CadShape:
    """Rescale ``shape`` by ``scale = (sx, sy, sz)`` (or a scalar) about its centroid.

    Preserves the shape's center of mass — scaling is performed about it.
    """
    require_cad()
    from .cad import transform_geometry  # noqa: PLC0415

    if isinstance(scale, (int, float)):
        sx = sy = sz = float(scale)
    else:
        sx, sy, sz = (float(s) for s in scale)

    center = shape.center()
    cx, cy, cz = center.x, center.y, center.z

    # translate(-c) → scale about origin → translate(+c), as a single 3x4 affine
    matrix = np.array(
        [
            [sx, 0.0, 0.0, cx - sx * cx],
            [0.0, sy, 0.0, cy - sy * cy],
            [0.0, 0.0, sz, cz - sz * cz],
        ],
        dtype=np.float64,
    )
    return transform_geometry(shape, matrix)


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

    :param phases: list of phases (must be CAD-backed)
    :param cut_obj: cutting object

    :return phase_cut: final result
    """
    require_cad()
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut  # noqa: PLC0415

    from .phase import Phase  # noqa: PLC0415

    phase_cut: list[Phase] = []

    for phase in phases:
        if phase.is_empty:
            continue
        cut = CadShape(BRepAlgoAPI_Cut(phase.cad.wrapped, cut_obj.wrapped).Shape())
        if len(cut.solids()) > 0:
            phase_cut.append(Phase.from_cad(cut))
    return phase_cut


def cut_phase_by_shape_list(phase_to_cut: Phase, shapes: list[CadShape]) -> Phase:
    """Cut a phase by a list of shapes.

    :param phase_to_cut: phase to cut (must be CAD-backed)
    :param shapes: list of cutting shapes

    :return resultCut: cut phase
    """
    require_cad()
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut  # noqa: PLC0415

    from .phase import Phase  # noqa: PLC0415

    if phase_to_cut.is_empty:
        err_msg = "phase_to_cut is empty"
        raise ValueError(err_msg)
    result = phase_to_cut.cad
    for shape in shapes:
        result = CadShape(BRepAlgoAPI_Cut(result.wrapped, shape.wrapped).Shape())
    return Phase.from_cad(result)


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
    """Cut list of phases in the given order (or reverse) and fuse them.

    :param phases: list of phases to cut (each must be CAD-backed)
    :param reverse_order: order for cutting shapes;
        when ``True``, the last shape of the list is not cut

    :return: list of phases
    """
    from .phase import Phase  # noqa: PLC0415

    shapes = [phase.cad for phase in phases]
    cutted_shapes = cut_shapes(shapes, reverse_order=reverse_order)

    return [Phase.from_cad(shape) for shape in cutted_shapes]


def _split_cad_in_grid(
    cad: CadShape,
    rve: Rve,
    grid: list[int],
) -> list[TopoDS_Shape]:
    """Split a CAD shape into solids per (grid-1)^3 interior cell planes."""
    require_cad()
    from .cad import enumerate_solids, make_plane_face, split_shape  # noqa: PLC0415

    result: list[TopoDS_Shape] = []
    for solid in enumerate_solids(cad):
        current = CadShape(solid)
        for axis in range(3):
            direction = tuple(int(axis == i) for i in range(3))
            coords = np.linspace(
                start=rve.min_point[axis],
                stop=rve.max_point[axis],
                num=grid[axis],
                endpoint=False,
            )[1:]
            for pos in coords:
                base_pnt = tuple(float(pos) * direction[k] for k in range(3))
                plane = make_plane_face(base_pnt, direction)
                current = split_shape(current, plane)
        result.extend(enumerate_solids(current))
    return result


def raster_phase(
    phase: Phase,
    rve: Rve,
    grid: list[int],
    *,
    phase_per_raster: bool = True,
) -> Phase | list[Phase]:
    """Raster a phase according to the RVE divided by the given grid.

    Each solid in the phase is split by the (grid-1) interior planes per
    axis (CAD-backed phases only).  When ``phase_per_raster`` is true,
    each non-empty grid cell becomes one :class:`Phase`; otherwise the
    split sub-solids are fused into one new :class:`Phase`.

    :param phase: CAD-backed phase to raster
    :param rve: RVE divided by the given grid
    :param grid: number of divisions in each direction ``[x, y, z]``
    :param phase_per_raster: if True, returns list of phases

    :return: Phase or list of Phases
    """
    require_cad()
    from OCP.BRepGProp import BRepGProp  # noqa: PLC0415
    from OCP.GProp import GProp_GProps  # noqa: PLC0415

    from .cad import make_compound_from_solids  # noqa: PLC0415
    from .phase import Phase  # noqa: PLC0415

    if phase.is_empty:
        err_msg = "Cannot raster an empty phase"
        raise ValueError(err_msg)

    solids = _split_cad_in_grid(phase.cad, rve, grid)

    if not phase_per_raster:
        return Phase.from_cad(make_compound_from_solids(solids))

    grid_arr = np.array(grid)
    buckets: list[list[TopoDS_Shape]] = [[] for _ in range(int(np.prod(grid_arr)))]
    for solid in solids:
        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(solid, props)
        com = props.CentreOfMass()
        center = np.array([com.X(), com.Y(), com.Z()])
        i, j, k = np.floor(
            grid_arr * (center - rve.min_point) / rve.dim,
        ).astype(int)
        ind = i + grid_arr[0] * j + grid_arr[0] * grid_arr[1] * k
        buckets[int(ind)].append(solid)
    return [
        Phase.from_cad(make_compound_from_solids(group)) for group in buckets if group
    ]


def repeat_shape(unit_geom: CadShape, rve: Rve, grid: tuple[int, int, int]) -> CadShape:
    """Repeat unit geometry in each direction according to the given grid.

    :param unit_geom: Shape to repeat
    :param rve: RVE of the geometry to repeat
    :param grid: list of number of geometry repetitions in each direction

    :return: CadShape of the repeated geometry
    """
    require_cad()
    from .cad import make_compound_from_solids, translate_solid  # noqa: PLC0415

    center = np.array(unit_geom.center().to_tuple())
    copies = []
    for idx in np.ndindex(*grid):
        pos = center - rve.dim * (0.5 * np.array(grid) - 0.5 - np.array(idx))
        copies.append(translate_solid(unit_geom.wrapped, pos))
    return make_compound_from_solids(copies)


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
