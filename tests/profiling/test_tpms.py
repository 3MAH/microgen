import microgen


def test_profiling_generate_gyroid() -> None:
    geometry = microgen.Tpms(
        surface_function=microgen.tpms.gyroid,
        type_part="sheet",
        thickness=0.1,
        repeat_cell=1,
    )
    geometry.generate(resolution=20)