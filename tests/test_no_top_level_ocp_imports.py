"""Lint-style test: no top-level OCP imports outside the CAD boundary.

The ``[cad]`` extra (``cadquery-ocp-novtk``) is strictly optional.  ``import
microgen`` must succeed without it, which means no module reachable from
``microgen.__init__`` may execute ``import OCP`` / ``from OCP …`` at
module-load time.  Every OCP import must be either:

  * inside a function/method body (lazy);
  * guarded by ``if TYPE_CHECKING:`` (only consumed by type checkers); or
  * located inside the CAD boundary (``microgen/cad.py`` today; the
    ``microgen/cad/`` subpackage after PR 4).

This test walks every ``.py`` under ``microgen/`` and AST-checks the
top level for stray OCP imports.  Failing here is what catches a future
refactor that accidentally drops an eager ``from OCP.… import …`` at
module scope and silently breaks ``pip install microgen`` (no extras).
"""
# ruff: noqa: S101

from __future__ import annotations

import ast
from pathlib import Path

MICROGEN_ROOT = Path(__file__).resolve().parent.parent / "microgen"

# Files allowed to have top-level OCP imports (the CAD boundary).  Paths are
# relative to ``microgen/``.  Extend this list as the CAD subpackage grows
# (PR 4 splits ``cad.py`` into a ``cad/`` subpackage).
_CAD_BOUNDARY = {
    "cad.py",
}


def _is_type_checking_block(node: ast.stmt) -> bool:
    """True if ``node`` is ``if TYPE_CHECKING:`` (any common spelling)."""
    if not isinstance(node, ast.If):
        return False
    test = node.test
    if isinstance(test, ast.Name) and test.id == "TYPE_CHECKING":
        return True
    if (
        isinstance(test, ast.Attribute)
        and isinstance(test.value, ast.Name)
        and test.value.id == "typing"
        and test.attr == "TYPE_CHECKING"
    ):
        return True
    return False


def _top_level_imports(tree: ast.Module) -> list[ast.stmt]:
    """Return module-level import statements, skipping ``if TYPE_CHECKING:``."""
    out: list[ast.stmt] = []
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            out.append(node)
        elif _is_type_checking_block(node):
            # Skip — these imports never execute at runtime.
            continue
    return out


def _imports_ocp(node: ast.stmt) -> bool:
    if isinstance(node, ast.ImportFrom):
        return (node.module or "").split(".", 1)[0] == "OCP"
    if isinstance(node, ast.Import):
        return any(alias.name.split(".", 1)[0] == "OCP" for alias in node.names)
    return False


def test_no_top_level_ocp_imports_outside_cad_boundary() -> None:
    """Walk microgen/*.py and assert no eager OCP import escapes the CAD boundary."""
    offenders: list[str] = []
    for py in MICROGEN_ROOT.rglob("*.py"):
        rel = py.relative_to(MICROGEN_ROOT).as_posix()
        if rel in _CAD_BOUNDARY:
            continue
        tree = ast.parse(py.read_text(encoding="utf-8"), filename=str(py))
        for node in _top_level_imports(tree):
            if _imports_ocp(node):
                offenders.append(f"{rel}:{node.lineno}")
    assert not offenders, (
        "Top-level OCP imports outside the CAD boundary "
        f"({sorted(_CAD_BOUNDARY)}): {offenders}. "
        "Move them inside a function body or under "
        "``if TYPE_CHECKING:``."
    )
