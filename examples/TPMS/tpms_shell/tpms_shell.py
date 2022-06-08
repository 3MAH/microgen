from microgen import Rve, Tpms, repeatGeometry, tpms, Phase
import cadquery as cq

shell = (
    cq.Workplane("front").box(3, 3, 3)
                         .faces("+Z or -X or +X")
                         .shell(0.1)
                         .translate((1., 1., 1.))
)
shell = shell.val()

rve = Rve(1, 1, 1)
geometry = Tpms(center=(0.5, 0.5, 0.5),
                surface_function=tpms.gyroid,
                type_part='sheet',
                thickness=0.2,
                path_data='data')
shape = geometry.generate()

inside = repeatGeometry(Phase(shape=shape), rve, {'x': 3, 'y': 3, 'z': 3})

compound = cq.Compound.makeCompound([shell, inside])
cq.exporters.export(compound, 'tpms_shell.stl')
