"""CAD backend install gate.

``require_cad()`` is the canonical way for lazy CAD code paths to verify
that the optional ``[cad]`` extra is available before touching OCP, and to
raise a user-friendly ``ImportError`` with install instructions otherwise.
"""

from __future__ import annotations

_INSTALL_HINT = (
    "microgen's CAD backend requires the OCP (OCCT) Python bindings. "
    "Install with:  pip install 'microgen[cad]'  "
    "(this pulls cadquery-ocp-novtk; on conda-forge use `ocp` instead)."
)


def require_cad() -> None:
    """Raise :class:`ImportError` if the CAD backend (OCP) is not importable."""
    try:
        import OCP  # noqa: F401, PLC0415
    except ImportError as err:
        raise ImportError(_INSTALL_HINT) from err
