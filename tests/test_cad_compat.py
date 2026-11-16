"""CAD wrapper contracts shared by OCP 7 and OCP 8."""

import pytest

pytest.importorskip("OCP")

from microgen.cad import make_box
from microgen.cad.topo import (
    enumerate_solids,
    make_plane_face,
    solid_center,
    split_shape,
)


def test_bounding_box_preserves_translated_extents():
    box = make_box((2, 4, 6), (3, -2, 5)).bounding_box()
    assert (box.xmin, box.ymin, box.zmin) == pytest.approx((2, -4, 2), abs=1e-6)
    assert (box.xmax, box.ymax, box.zmax) == pytest.approx((4, 0, 8), abs=1e-6)


@pytest.mark.parametrize("multiple_tools", [False, True])
def test_split_shape_retains_all_arguments_and_tools(multiple_tools):
    box = make_box((3, 2, 2), (0, 0, 0))
    tools = [make_plane_face((-0.5, 0, 0), (1, 0, 0))]
    if multiple_tools:
        tools.append(make_plane_face((0.5, 0, 0), (1, 0, 0)))
    try:
        result = split_shape(box, tools if multiple_tools else tools[0])
    except ImportError as error:
        pytest.fail(f"Installed OCP cannot split CAD shapes: {error}")
    solids = enumerate_solids(result)
    assert len(solids) == len(tools) + 1
    assert result.volume() == pytest.approx(12)
    centers = sorted(solid_center(solid)[0] for solid in solids)
    expected = [-1, 0, 1] if multiple_tools else [-1, 0.5]
    assert centers == pytest.approx(expected)


@pytest.mark.parametrize("cell_size", [1.0, 2.0])
def test_thin_lattice_fusion_preserves_small_solid_volume(cell_size):
    from microgen.shape.strut_lattice import Diamond

    lattice = Diamond(strut_radius=0.001, cell_size=cell_size)
    assert 0 < lattice.volume < 0.001 * cell_size**3
