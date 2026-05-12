"""Top-level pytest configuration.

Tests that hit ``microgen.cad.require_cad()`` raise ``ImportError`` when
the optional ``[cad]`` extra is not installed. CI's ``cad`` job has OCP;
local installs without ``[cad]`` don't. We catch that specific
``ImportError`` here and convert it into a pytest skip, so the dev-env
test run reports clean skips instead of dozens of identical noise
failures. Tests that don't touch CAD are unaffected.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from collections.abc import Generator


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item: pytest.Item) -> Generator[None, None, None]:
    """Convert OCP-missing ``ImportError`` / ``ModuleNotFoundError`` into a skip.

    Accepts two shapes:
    - ``ModuleNotFoundError("No module named 'OCP'")`` from internal
      lazy imports of OCP inside microgen.
    - ``ImportError("microgen's CAD backend requires …")`` raised by
      :func:`microgen.cad.require_cad`.
    """
    outcome = yield
    excinfo = outcome.excinfo if hasattr(outcome, "excinfo") else None
    if excinfo is None:
        return
    exc_type, exc_value, _tb = excinfo
    if not issubclass(exc_type, ImportError):
        return
    message = str(exc_value)
    name = getattr(exc_value, "name", "") or ""
    if name == "OCP" or "OCP" in message or "microgen's CAD backend" in message:
        pytest.skip(f"OCP (CAD backend) not installed: {exc_value}")
