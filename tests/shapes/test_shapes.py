from typing import Type

import numpy as np
import pytest

import microgen
from microgen.shape import (
    BasicGeometry,
    Box,
    Capsule,
    Cylinder,
    Ellipsoid,
    ExtrudedPolygon,
    Polyhedron,
    Sphere,
)


def test_shapes():
    rve = microgen.Rve(dim=1)

    elem = microgen.newGeometry(
        shape="Ellipsoid", param_geom={"a_x": 0.15, "a_y": 0.31, "a_z": 0.4}
    )
    ellipsoid = elem.generate()
    phase = microgen.Phase(shape=ellipsoid)
    phase.centerOfMass
    phase.centerOfMass
    phase.getCenterOfMass(compute=False)
    phase.inertiaMatrix
    phase.inertiaMatrix
    phase.getInertiaMatrix(compute=False)
    phase.solids
    phase.shape
    phase.translate((1, 0, 0))
    phase.translate(np.array([0, 1, 1]))
    phase.rescale(1.5)
    phase.repeat(rve, (1, 2, 1))
    phase.rasterize(rve, [2, 2, 2], phasePerRaster=True)
    phase.rasterize(rve, [2, 2, 2], phasePerRaster=False)

    void_phase = microgen.Phase()
    void_phase.shape
    void_phase.solids

    elem = microgen.newGeometry(shape="Sphere", param_geom={"radius": 0.15})
    elem = microgen.newGeometry(
        shape="Box", param_geom={"dim_x": 0.15, "dim_y": 0.31, "dim_z": 0.4}
    )
    elem = microgen.newGeometry(
        shape="Capsule", param_geom={"height": 0.5, "radius": 0.1}
    )
    elem = microgen.newGeometry(
        shape="Cylinder", param_geom={"height": 0.5, "radius": 0.1}
    )
    elem = microgen.newGeometry(
        shape="ExtrudedPolygon",
        param_geom={"listCorners": [(0, 0), (0, 1), (1, 1), (1, 0)], "height": 0.3},
    )
    dic = microgen.shape.polyhedron.read_obj(
        "examples/BasicShapes/platon/tetrahedron.obj"
    )
    microgen.newGeometry(shape="Polyhedron", param_geom={"dic": dic})

    with pytest.raises(ValueError):
        microgen.newGeometry(shape="fake", param_geom={"fake": 0})


@pytest.mark.parametrize(
    "shape",
    [Box, Capsule, Cylinder, Ellipsoid, ExtrudedPolygon, Polyhedron, Sphere],
)
def test_shape_cad_and_vtk_volume_must_correspond(shape: Type[BasicGeometry]):
    geom = shape()
    shape_cad = geom.generate()
    shape_vtk = geom.generateVtk()

    volume_cad = shape_cad.Volume()

    assert volume_cad > 0
    assert np.isclose(volume_cad, shape_vtk.volume, rtol=1e-2)
