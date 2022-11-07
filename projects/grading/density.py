from microgen import Tpms, tpms, fuseShapes

import numpy as np
import pyvista as pv
import cadquery as cq

import itertools

def createSurfaces(
    geometry,
    isovalues: list[float],
    resolution: int = 20,
    smoothing: int = 100,
) -> list[cq.Shell]:
    """
    Create TPMS surfaces for the corresponding isovalue, return a list of cq.Shell

    :param numsber_surfaces: number of surfaces
    :param isovalues: height isovalues of the given tpms function
    :param resolution: surface file name
    :param smoothing: smoothing loop iterations
    """

    grid, surface = geometry.computeSurfaceFunction(geometry.surface_function, resolution)

    shells = []
    for isovalue in isovalues:
        mesh = geometry.computeSurfaceVtk(grid, surface, isovalue, smoothing)
        shells.append(geometry.createShell(mesh))
    return shells

def generate(
    geometry,
    resolution: int = 20,
    smoothing: int = 100,
) -> cq.Shape:
    """
    Creates thick TPMS geometry (sheet or skeletal part) from surface

    :param resolution: surface file name
    :param smoothing: smoothing loop iterations
    """

    isovalues = [
        -geometry.thickness / 2.,
        -geometry.thickness / 6.0,
        geometry.thickness / 6.0,
        geometry.thickness / 2.,
    ]
    shells = createSurfaces(geometry, isovalues=isovalues, resolution=resolution, smoothing=smoothing)

    face_cut_tp = shells[2]
    face_cut_tm = shells[1]
    face_cut_p = shells[3]
    face_cut_m = shells[0]

    box_wp = cq.Workplane("front").box(1, 1, 1)

    boxCut_wp = box_wp.split(face_cut_p)
    boxCut_wp = boxCut_wp.split(face_cut_m)

    boxWorkplanes: list[cq.Workplane] = boxCut_wp.solids().all()

    listShapes: list[tuple[int, cq.Shape]] = []
    for wp in boxWorkplanes:
        n_shapes = wp.split(face_cut_tp).split(face_cut_tm).solids().size()
        shape = wp.val()
        listShapes.append((n_shapes, shape))

    if geometry.type_part == "sheet":
        sheet = [shape for (n_shapes, shape) in listShapes if n_shapes > 1]
        to_fuse = [cq.Shape(shape.wrapped) for shape in sheet]
        shape = fuseShapes(to_fuse, True)
    elif geometry.type_part == "skeletal":
        skeletal = [shape for (n_shapes, shape) in listShapes if n_shapes == 1]
        shape = skeletal[0]
        # to_fuse = [cq.Shape(shape.wrapped) for shape in skeletal]
        # shape = fuseShapes(to_fuse, False)
    return shape


geometry = Tpms(
    surface_function=tpms.gyroid,
    type_part="sheet",
    thickness=0.2,
    repeat_cell=1,
    cell_size=1
)

shape = generate(geometry)
shape = pv.PolyData(shape.toVtkPolyData(tolerance=0.01, angularTolerance=0.1, normals=True))
p = pv.Plotter()
p.add_mesh(shape, color='w')
p.camera_position = 'xy'
p.parallel_projection = True 
p.show_axes()
p.show()