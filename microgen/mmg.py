import os

def mmg2d(d=None, h=None, m=None, v=None, val=None, default=None,        \
          input=None, output=None, solution=None, metric=None,           \
          A=None, ar=None, hausd=None, hgrad=None, hmax=None,            \
          hmin=None, hsiz=None, lag=None, ls=None, _3dMedit=None,        \
          noinsert=None, nomove=None, nosurf=None, noswap=None, nr=None, \
          nreg=None, nsd=None, optim=None, opnbdy=None, rmc=None):
    pass

def mmgs(d=None, h=None, m=None, v=None, val=None, default=None,        \
          input=None, output=None, solution=None, metric=None,          \
          A=None, ar=None, hausd=None, hgrad=None, hmax=None, hmin=None,\
          hsiz=None, ls=None, noinsert=None, nomove=None, nosurf=None,  \
          noswap=None, nr=None, nreg=None, nsd=None, optim=None, rn=None):
    pass

def mmg3d(d=None, h=None, m=None, v=None, val=None, default=None,\
          input=None, output=None, solution=None, metric=None,   \
          A=None, ar=None, octree=None, hausd=None, hgrad=None,  \
          hmax=None, hmin=None, hsiz=None, lag=None, ls=None,    \
          nofem=None, noinsert=None, nomove=None, nosurf=None,   \
          noswap=None, nr=None, nreg=None, nsd=None, optim=None, \
          optimLES=None, opnbdy=None, rmc=None, rn=None):
    cmd = "mmg3d_O3"
    if d:
        cmd += " -d "
    if h:
        cmd += " -h "
    if m:
        if m == True:
            m = ""
        cmd += " -m " + str(m)
    if v:
        if v == True:
            v = 1
        cmd += " -v " + str(v)
    if val:
        cmd += " -val "
    if default:
        cmd += " -default "
    if input:
        cmd += " -in " + input
    if output:
        cmd += " -out " + output
    if solution:
        cmd += " -sol " + solution
    if metric:
        cmd += " -met " + metric
    if A:
        cmd += " -A "
    if ar:
        cmd += " -ar " + str(ar)
    if octree:
        cmd += " -octree " + str(octree)
    if hausd:
        cmd += " -hausd " + str(hausd)
    if hgrad:
        cmd += " -hgrad " + str(hgrad)
    if hmax:
        cmd += " -hmax " + str(hmax)
    if hmin:
        cmd += " -hmin " + str(hmin)
    if hsiz:
        cmd += " -hsiz " + str(hsiz)
    if lag or lag == 0:
        if lag == True:
            lag = 0
        cmd += " -lag " + str(lag)
    if ls or ls == 0:
        ls_value = ls
        if ls_value == True:
            ls_value = 0
        cmd += " -ls " + str(ls_value)
    if nofem:
        cmd += " -nofem "
    if noinsert:
        cmd += " -noinsert "
    if nomove:
        cmd += " -nomove "
    if nosurf:
        cmd += " -nosurf "
    if noswap:
        cmd += " -noswap "
    if nr:
        cmd += " -nr "
    if nreg:
        cmd += " -nreg " + str(nreg)
    if nsd:
        cmd += " -nsd " + str(nsd)
    if optim:
        cmd += " -optim "
    if optimLES:
        cmd += " -optimLES "
    if opnbdy:
        cmd += " -opnbdy "
    if rmc:
        cmd += " -rmc " + str(rmc)
    if rn:
        cmd += " -rn "

    print(cmd)
    os.system(cmd)