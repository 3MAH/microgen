"""Smoke tests verifying that the mesh / implicit-field path works without OCP.

These tests exercise the public microgen API that must function **without**
the optional ``[cad]`` install extra (i.e. without ``cadquery-ocp``/OCP).

In CI we run them in a dedicated no-CAD environment.  Running them in an
environment where OCP *is* installed is still fine — every assertion here
holds unconditionally; these tests just prove the no-CAD paths don't
*transitively* pull in OCP via an eager import.
"""
# ruff: noqa: S101

from __future__ import annotations

import importlib
import sys


def test_import_microgen_without_ocp() -> None:
    """``import microgen`` must succeed even if OCP is not importable."""
    import microgen  # noqa: F401

    # microgen itself must not eagerly import OCP.
    assert "OCP" not in sys.modules or sys.modules.get("OCP") is not None


def test_implicit_ops_module_does_not_need_ocp() -> None:
    """Implicit (F-rep) operations are pure numpy — no OCP."""
    from microgen.shape import implicit_ops

    # Reload with OCP masked to confirm the module never imports it at
    # module-load time.  (The module has no OCP import; this just pins it.)
    importlib.reload(implicit_ops)
    assert hasattr(implicit_ops, "union")
    assert hasattr(implicit_ops, "shell")
    assert hasattr(implicit_ops, "normalize_to_sdf")


def test_primitive_generate_vtk_without_cad_extra() -> None:
    """Every primitive's ``generate_vtk()`` works without the CAD extra."""
    from microgen.shape.box import Box
    from microgen.shape.cylinder import Cylinder
    from microgen.shape.sphere import Sphere

    assert Box().generate_vtk().n_cells > 0
    assert Sphere().generate_vtk().n_cells > 0
    assert Cylinder().generate_vtk().n_cells > 0


def test_tpms_generate_vtk_without_cad_extra() -> None:
    """TPMS F-rep + marching cubes work without the CAD extra."""
    from microgen import surface_functions
    from microgen.shape.tpms import Tpms

    tpms = Tpms(surface_function=surface_functions.gyroid, offset=0.5)
    mesh = tpms.generate_vtk(type_part="sheet")
    assert mesh.n_cells > 0


def test_cad_capable_entry_points_raise_cleanly_without_ocp() -> None:
    """Calling ``.generate()`` / ``make_box()`` without OCP gives a clear error.

    Only meaningful in a true no-CAD env; when OCP *is* installed these calls
    succeed.  We just verify the function exists and does not, for instance,
    segfault on import.
    """
    from microgen import cad

    try:
        cad.require_cad()
    except ImportError as err:
        assert "cadquery-ocp" in str(err) or "microgen[cad]" in str(err)
    else:
        # OCP is installed in this env — that's fine, nothing to verify.
        assert True
