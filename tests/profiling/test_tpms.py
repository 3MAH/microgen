import microgen


def test_profiling_generate_gyroid() -> None:
    geometry = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=0.3,
        repeat_cell=1,
    )
    geometry.generateVtk(type_part="sheet")