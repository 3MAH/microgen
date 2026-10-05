"""Transform contract of :class:`Phase`: pivot, ``iso``, ``period``, cached grid, CAD."""

import numpy as np
import pytest
import pyvista as pv
from scipy.spatial.transform import Rotation

from microgen import Box, Phase, Rve, Sphere, Tpms, surface_functions

# ruff: noqa: S101

_ROT_Z90 = Rotation.from_euler("z", 90, degrees=True)
_ROT_GENERIC = Rotation.from_euler("xyz", (30, 20, 10), degrees=True)
_PTS = np.random.default_rng(0).uniform(-2.0, 2.0, size=(3, 500))


def _sphere_phase(**kwargs) -> Phase:
    return Phase.from_shape(Sphere(radius=0.5, center=(0.3, 0.0, 0.0)), **kwargs)


@pytest.mark.parametrize(
    ("verb", "args"),
    [
        ("translate", ((1.0, 0.0, -0.5),)),
        ("rotate", (_ROT_GENERIC,)),
        ("scale", ((2.0, 1.0, 0.5),)),
    ],
)
def test_inplace_matches_copy(verb, args) -> None:
    ref = getattr(_sphere_phase(), verb)(*args)
    phase = _sphere_phase()
    out = getattr(phase, verb)(*args, inplace=True)
    assert out is phase
    np.testing.assert_allclose(out.field(*_PTS), ref.field(*_PTS))
    assert out.bounds == ref.bounds
    assert out.iso == ref.iso


def test_scale_defaults_to_origin_pivot() -> None:
    phase = _sphere_phase()
    out = phase.scale(2.0)
    # The center (0.3, 0, 0) moves to (0.6, 0, 0); the radius doubles.
    assert out.field(np.array([0.6]), np.array([0.0]), np.array([0.0]))[0] == (
        pytest.approx(-1.0)
    )


def test_scale_about_point_keeps_point_fixed() -> None:
    phase = _sphere_phase()
    p = np.array([[0.3], [0.0], [0.0]])
    out = phase.scale(3.0, point=p.ravel())
    np.testing.assert_allclose(out.field(*p), 3.0 * phase.field(*p))


def test_scale_rescales_iso() -> None:
    """With ``iso != 0`` the solid scales exactly: volume ratio = det(S)."""
    phase = Phase.from_shape(Sphere(radius=0.4), iso=0.1, resolution=60)
    out = phase.scale((2.0, 2.0, 2.0))
    assert out.iso == pytest.approx(0.2)
    # {f < iso} of the scaled phase is the sphere of radius 2 * 0.5.
    assert out.field(np.array([0.99]), np.array([0.0]), np.array([0.0]))[0] < out.iso
    assert out.field(np.array([1.01]), np.array([0.0]), np.array([0.0]))[0] > out.iso


def test_scale_rejects_non_positive() -> None:
    with pytest.raises(ValueError, match="strictly positive"):
        _sphere_phase().scale(-2.0)


def test_period_rules() -> None:
    tpms = Tpms(
        surface_function=surface_functions.gyroid,
        offset=0.3,
        cell_size=(1.0, 2.0, 1.0),
    )
    phase = Phase.from_shape(tpms)
    assert phase.translate((0.5, 0.0, 0.0)).period == (1.0, 2.0, 1.0)
    assert phase.scale((2.0, 1.0, 3.0)).period == (2.0, 2.0, 3.0)
    assert phase.rotate(_ROT_Z90).period == (2.0, 1.0, 1.0)
    assert phase.rotate(_ROT_GENERIC).period is None


def _non_cubic_grid() -> pv.StructuredGrid:
    x, y, z = np.meshgrid(
        np.linspace(-1.0, 1.0, 21),
        np.linspace(-1.0, 1.0, 31),
        np.linspace(-0.5, 0.5, 11),
        indexing="ij",
    )
    sg = pv.StructuredGrid(x, y, z)
    sg["implicit"] = (np.sqrt(x**2 + y**2 + z**2) - 0.43).ravel(order="F")
    return sg


def test_from_grid_non_cubic_moments_and_scale() -> None:
    phase = Phase.from_grid(_non_cubic_grid())
    com = phase.center_of_mass
    np.testing.assert_allclose(com, 0.0, atol=1e-12)
    out = phase.scale(2.0)
    assert out.grid() is not None
    np.testing.assert_allclose(out.grid().bounds, 2.0 * np.array(phase.grid().bounds))
    # Same number of inside nodes, cell volume x 8.
    ratio = out.inertia_matrix[0, 0] / phase.inertia_matrix[0, 0]
    assert ratio == pytest.approx(2.0**5)


def test_rotated_from_grid_keeps_inertia_trace() -> None:
    """The trace of the inertia tensor about the origin is rotation-invariant."""
    phase = Phase.from_grid(_non_cubic_grid())
    rotated = phase.rotate(_ROT_GENERIC)
    assert np.trace(rotated.inertia_matrix) == pytest.approx(
        np.trace(phase.inertia_matrix), rel=1e-10
    )
    assert len(rotated.pieces) == len(phase.pieces) == 1
    assert rotated.pieces[0].volume == pytest.approx(phase.pieces[0].volume)


def test_from_grid_custom_scalar_name() -> None:
    sg = _non_cubic_grid()
    sg["sdf"] = sg["implicit"]
    del sg.point_data["implicit"]
    phase = Phase.from_grid(sg, scalars="sdf")
    np.testing.assert_allclose(phase.center_of_mass, 0.0, atol=1e-12)


def test_cad_uniform_scale_keeps_analytic_surface() -> None:
    pytest.importorskip("OCP")
    from OCP.BRep import BRep_Tool
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    phase = Phase.from_cad(Sphere(radius=0.5, center=(1.0, 0.0, 0.0)).generate_cad())
    out = phase.scale(2.0)
    face = TopExp_Explorer(out.cad.wrapped, TopAbs_FACE).Current()
    surface = BRep_Tool.Surface_s(TopoDS.Face_s(face))
    assert surface.DynamicType().Name() == "Geom_SphericalSurface"
    assert out.cad.volume() == pytest.approx(8.0 * phase.cad.volume(), rel=1e-9)
    np.testing.assert_allclose(out.center_of_mass, (2.0, 0.0, 0.0), atol=1e-9)


def test_cad_anisotropic_scale_volume() -> None:
    pytest.importorskip("OCP")
    phase = Phase.from_cad(Box(dim=(1.0, 1.0, 1.0)).generate_cad())
    out = phase.scale((1.0, 2.0, 3.0))
    assert out.cad.volume() == pytest.approx(6.0, rel=1e-9)


def test_tile_single_copy_stays_in_place() -> None:
    pytest.importorskip("OCP")
    phase = Phase.from_cad(
        Box(dim=(0.5, 0.5, 0.5), center=(0.2, 0.1, 0.0)).generate_cad()
    )
    tiled = phase.tile(Rve(dim=1.0), (1, 1, 1))
    np.testing.assert_allclose(tiled.center_of_mass, (0.2, 0.1, 0.0), atol=1e-9)


def test_generic_rotation_clips_to_rotated_domain() -> None:
    """A periodic field rotated by 45 degrees brings no material from outside."""
    tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.3)
    phase = Phase.from_shape(tpms)
    rotated = phase.rotate(Rotation.from_euler("z", 45, degrees=True))
    # Inside the new AABB but outside the rotated cell [-0.5, 0.5]^3.
    corner = (np.array([0.6]), np.array([0.3]), np.array([0.0]))
    assert rotated.field(*corner)[0] >= rotated.iso
    # Inside the rotated cell the field is the rotated field.
    p = np.array([[0.1], [0.05], [0.2]])
    back = Rotation.from_euler("z", -45, degrees=True).apply(p.T).T
    assert rotated.field(*p)[0] == pytest.approx(phase.field(*back)[0])


def test_inertia_is_about_center_of_mass() -> None:
    phase = Phase.from_shape(Sphere(radius=0.5), resolution=40)
    moved = phase.translate((1.0, 0.5, 0.0))
    np.testing.assert_allclose(moved.inertia_matrix, phase.inertia_matrix, atol=1e-12)


def test_inertia_field_matches_cad() -> None:
    pytest.importorskip("OCP")
    sphere = Sphere(radius=0.5, center=(1.0, 0.0, 0.0))
    field_phase = Phase.from_shape(sphere, resolution=80)
    cad_phase = Phase.from_cad(sphere.generate_cad())
    np.testing.assert_allclose(
        np.diag(field_phase.inertia_matrix),
        np.diag(cad_phase.inertia_matrix),
        rtol=0.05,
    )


def test_from_grid_rotated_grid_nearest_field() -> None:
    sg = _non_cubic_grid().rotate(
        Rotation.from_euler("z", 30, degrees=True), inplace=False
    )
    phase = Phase.from_grid(sg)
    pts = np.asarray(sg.points)
    np.testing.assert_array_equal(phase.field(*pts.T), np.asarray(sg["implicit"]))


def test_transform_leaves_original_cad_untouched() -> None:
    pytest.importorskip("OCP")
    phase = Phase.from_cad(Sphere(radius=0.5).generate_cad())
    phase.translate((2.0, 0.0, 0.0))
    np.testing.assert_allclose(
        phase.cad.center().to_tuple(), (0.0, 0.0, 0.0), atol=1e-9
    )


def test_mesh_backed_scale_and_rotate() -> None:
    phase = Phase.from_mesh(
        pv.Sphere(radius=0.5, theta_resolution=40, phi_resolution=40)
    )
    scaled = phase.scale(2.0)
    assert abs(scaled.surface_mesh().volume) == pytest.approx(
        8.0 * abs(phase.surface_mesh().volume), rel=1e-9
    )
    rotated = phase.rotate(_ROT_GENERIC, point=(1.0, 0.0, 0.0))
    np.testing.assert_allclose(
        rotated.surface_mesh().center,
        _ROT_GENERIC.apply((-1.0, 0.0, 0.0)) + (1.0, 0.0, 0.0),
        atol=1e-6,
    )


def test_cad_reflection_keeps_positive_volume() -> None:
    pytest.importorskip("OCP")
    from microgen.cad import transform_geometry

    box = Box(dim=(1.0, 2.0, 3.0)).generate_cad()
    mirrored = transform_geometry(
        box, np.hstack([np.diag([-1.0, 1.0, 1.0]), np.zeros((3, 1))])
    )
    assert mirrored.volume() == pytest.approx(6.0, rel=1e-9)
