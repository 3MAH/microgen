import os

def mmg3d(mesh, output, ls, hsiz):
    cmd = "mmg3d_O3 "
    cmd += mesh
    if output:
        cmd += " -out " + output
    if ls:
        cmd += " -ls "
    if hsiz:
        cmd += " -hsiz " + str(hsiz)

    print(cmd)
    os.system(cmd)

def mmg2d():
    pass