import subprocess

def mmg2d(d=None, h=None, m=None, v=None, val=None, default=None,        \
          input=None, output=None, solution=None, metric=None,           \
          A=None, ar=None, hausd=None, hgrad=None, hmax=None,            \
          hmin=None, hsiz=None, lag=None, ls=None, _3dMedit=None,        \
          noinsert=None, nomove=None, nosurf=None, noswap=None, nr=None, \
          nreg=None, nsd=None, optim=None, opnbdy=None, rmc=None):
    cmd = ["mmg2d_O3"]
    if d:
        cmd.append("-d")
    if h:
        cmd.append("-h")
    if m:
        if m == True:
            m = ""
        cmd.append("-m " + str(m))
    if v:
        if v == True:
            v = 1
        cmd.append("-v " + str(v))
    if val:
        cmd.append("-val")
    if default:
        cmd.append("-default")
    if input:
        cmd.append("-in " + input)
    if output:
        cmd.append("-out " + output)
    if solution:
        cmd.append("-sol " + solution)
    if metric:
        cmd.append("-met " + metric)
    if A:
        cmd.append("-A")
    if ar:
        cmd.append("-ar " + str(ar))
    if hausd:
        cmd.append("-hausd " + str(hausd))
    if hgrad:
        cmd.append("-hgrad " + str(hgrad))
    if hmax:
        cmd.append("-hmax " + str(hmax))
    if hmin:
        cmd.append("-hmin " + str(hmin))
    if hsiz:
        cmd.append("-hsiz " + str(hsiz))
    if lag or lag == 0:
        if lag == True:
            lag = 0
        cmd.append("-lag " + str(lag))
    if ls or ls == 0:
        ls_value = ls
        if ls_value == True:
            ls_value = 0
        cmd.append("-ls " + str(ls_value))
    if _3dMedit:
        cmd.append("-3dMedit " + str(_3dMedit))
        
    if noinsert:
        cmd.append("-noinsert")
    if nomove:
        cmd.append("-nomove")
    if nosurf:
        cmd.append("-nosurf")
    if noswap:
        cmd.append("-noswap")
    if nr:
        cmd.append("-nr")
    if nreg:
        cmd.append("-nreg " + str(nreg))
    if nsd:
        cmd.append("-nsd " + str(nsd))
    if optim:
        cmd.append("-optim")
    if opnbdy:
        cmd.append("-opnbdy")
    if rmc:
        cmd.append("-rmc " + str(rmc))

    subprocess.run(cmd)

def mmgs(d=None, h=None, m=None, v=None, val=None, default=None,        \
          input=None, output=None, solution=None, metric=None,          \
          A=None, ar=None, hausd=None, hgrad=None, hmax=None, hmin=None,\
          hsiz=None, ls=None, noinsert=None, nomove=None, nosurf=None,  \
          noswap=None, nr=None, nreg=None, nsd=None, optim=None, rn=None):
    cmd = "mmgs_O3"
    if d:
        cmd.append("-d")
    if h:
        cmd.append("-h")
    if m:
        if m == True:
            m = ""
        cmd.append("-m " + str(m))
    if v:
        if v == True:
            v = 1
        cmd.append("-v " + str(v))
    if val:
        cmd.append("-val")
    if default:
        cmd.append("-default")
    if input:
        cmd.append("-in " + input)
    if output:
        cmd.append("-out " + output)
    if solution:
        cmd.append("-sol " + solution)
    if metric:
        cmd.append("-met " + metric)
    if A:
        cmd.append("-A")
    if ar:
        cmd.append("-ar " + str(ar))
    if hausd:
        cmd.append("-hausd " + str(hausd))
    if hgrad:
        cmd.append("-hgrad " + str(hgrad))
    if hmax:
        cmd.append("-hmax " + str(hmax))
    if hmin:
        cmd.append("-hmin " + str(hmin))
    if hsiz:
        cmd.append("-hsiz " + str(hsiz))
    if ls or ls == 0:
        ls_value = ls
        if ls_value == True:
            ls_value = 0
        cmd.append("-ls " + str(ls_value))
    if noinsert:
        cmd.append("-noinsert")
    if nomove:
        cmd.append("-nomove")
    if nosurf:
        cmd.append("-nosurf")
    if noswap:
        cmd.append("-noswap")
    if nr:
        cmd.append("-nr")
    if nreg:
        cmd.append("-nreg " + str(nreg))
    if nsd:
        cmd.append("-nsd " + str(nsd))
    if optim:
        cmd.append("-optim")
    if rn:
        cmd.append("-rn")

    print(cmd)
    subprocess.run(cmd)

def mmg3d(d=None, h=None, m=None, v=None, val=None, default=None,\
          input=None, output=None, solution=None, metric=None,   \
          A=None, ar=None, octree=None, hausd=None, hgrad=None,  \
          hmax=None, hmin=None, hsiz=None, lag=None, ls=None,    \
          nofem=None, noinsert=None, nomove=None, nosurf=None,   \
          noswap=None, nr=None, nreg=None, nsd=None, optim=None, \
          optimLES=None, opnbdy=None, rmc=None, rn=None):
    cmd = "mmg3d_O3"
    if d:
        cmd.append("-d")
    if h:
        cmd.append("-h")
    if m:
        if m == True:
            m = ""
        cmd.append("-m " + str(m))
    if v:
        if v == True:
            v = 1
        cmd.append("-v " + str(v))
    if val:
        cmd.append("-val")
    if default:
        cmd.append("-default")
    if input:
        cmd.append("-in " + input)
    if output:
        cmd.append("-out " + output)
    if solution:
        cmd.append("-sol " + solution)
    if metric:
        cmd.append("-met " + metric)
    if A:
        cmd.append("-A")
    if ar:
        cmd.append("-ar " + str(ar))
    if octree:
        cmd.append("-octree " + str(octree))
    if hausd:
        cmd.append("-hausd " + str(hausd))
    if hgrad:
        cmd.append("-hgrad " + str(hgrad))
    if hmax:
        cmd.append("-hmax " + str(hmax))
    if hmin:
        cmd.append("-hmin " + str(hmin))
    if hsiz:
        cmd.append("-hsiz " + str(hsiz))
    if lag or lag == 0:
        if lag == True:
            lag = 0
        cmd.append("-lag " + str(lag))
    if ls or ls == 0:
        ls_value = ls
        if ls_value == True:
            ls_value = 0
        cmd.append("-ls " + str(ls_value))
    if nofem:
        cmd.append("-nofem")
    if noinsert:
        cmd.append("-noinsert")
    if nomove:
        cmd.append("-nomove")
    if nosurf:
        cmd.append("-nosurf")
    if noswap:
        cmd.append("-noswap")
    if nr:
        cmd.append("-nr")
    if nreg:
        cmd.append("-nreg " + str(nreg))
    if nsd:
        cmd.append("-nsd " + str(nsd))
    if optim:
        cmd.append("-optim")
    if optimLES:
        cmd.append("-optimLES")
    if opnbdy:
        cmd.append("-opnbdy")
    if rmc:
        cmd.append("-rmc " + str(rmc))
    if rn:
        cmd.append("-rn")

    subprocess.run(cmd)