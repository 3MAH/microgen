import microgen
import pytest
from microgen.external import MmgError


def test_run_mmg2d_with_no_arguments_must_raise_MmgError() -> None:
    with pytest.raises(MmgError):
        microgen.external.Mmg.mmg2d()


def test_run_mmgs_with_no_arguments_must_raise_MmgError() -> None:
    with pytest.raises(MmgError):
        microgen.external.Mmg.mmgs()


def test_run_mmg3d_with_no_arguments_must_raise_MmgError() -> None:
    with pytest.raises(MmgError):
        microgen.external.Mmg.mmg3d()
