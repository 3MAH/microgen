import microgen
import cadquery as cq

print("generate shell")
shell = (cq.Workplane("front")
            .box(3, 3, 3)
            .faces("+Z or -X or +X")
            .shell(0.1)
            .translate((1., 1., 1.))
          )
shell = shell.val()
# cq.exporters.export(shell, 'shell.step')
# cq.exporters.export(shell, 'shell.stl')

print("generate gyroid")
rve = microgen.Rve(1, 1, 1)
phase = microgen.BasicGeometry(number=0, shape='tpms',
                               xc=0.5, yc=0.5, zc=0.5,
                               psi=0., theta=0., phi=0.,
                               param_geom={"surface_function": microgen.shape.tpms.gyroid,
                                           "type_part": 'sheet',
                                           "thickness": 0.2},
                               path_data='data')
phase.geometry.createSurfaces(rve=rve,
                              sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03,
                              path_data='data')
gyroid = phase.generate(rve=rve)

print("repeat geometry")
inside = microgen.repeatGeometry(gyroid, rve, {'x': 3, 'y': 3, 'z': 3})
# cq.exporters.export(inside, 'inside.step')
# cq.exporters.export(inside, 'inside.stl')

compound = cq.Compound.makeCompound([shell, inside])
# cq.exporters.export(compound, 'tpms_shell.step')
cq.exporters.export(compound, 'tpms_shell.stl')