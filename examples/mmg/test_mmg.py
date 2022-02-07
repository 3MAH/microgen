import microgen

mesh = "initialmesh.mesh"
microgen.mmg.mmg3d(input=mesh, output="finalmesh.mesh", ls=True, hsiz=0.03)
