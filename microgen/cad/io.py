"""STEP file I/O (read).

STEP/BREP/STL *export* is exposed as methods on :class:`CadShape` itself
(``export_step``, ``export_brep``, ``export_stl``) — see ``cad/shape.py``.
This module holds the standalone ``import_step`` reader.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ._install import require_cad
from .shape import CadShape

if TYPE_CHECKING:
    from pathlib import Path


def import_step(path: str | Path) -> CadShape:
    """Import a STEP file and return the resulting :class:`CadShape`.

    Multi-root STEP files are merged into a single ``TopoDS_Compound``.
    """
    require_cad()
    from OCP.IFSelect import IFSelect_RetDone  # noqa: PLC0415
    from OCP.STEPControl import STEPControl_Reader  # noqa: PLC0415

    reader = STEPControl_Reader()
    status = reader.ReadFile(str(path))
    if status != IFSelect_RetDone:
        err_msg = f"STEP read failed for {path!r} with status {status!r}"
        raise RuntimeError(err_msg)
    reader.TransferRoots()
    return CadShape(reader.OneShape())
