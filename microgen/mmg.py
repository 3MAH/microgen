import subprocess


def mmg2d(
    d=None,
    h=None,
    m=None,
    v=None,
    val=None,
    default=None,
    input=None,
    output=None,
    solution=None,
    metric=None,
    A=None,
    ar=None,
    hausd=None,
    hgrad=None,
    hmax=None,
    hmin=None,
    hsiz=None,
    lag=None,
    ls=None,
    _3dMedit=None,
    noinsert=None,
    nomove=None,
    nosurf=None,
    noswap=None,
    nr=None,
    nreg=None,
    nsd=None,
    optim=None,
    opnbdy=None,
    rmc=None,
):
    cmd = ["mmg2d_O3"]
    if d:
        cmd.append("-d")
    if h:
        cmd.append("-h")
    if m:
        if isinstance(m, bool):
            m = ""
        cmd.append("-m")
        cmd.append(str(m))
    if v:
        if isinstance(v, bool):
            v = 1
        cmd.append("-v")
        cmd.append(str(v))
    if val:
        cmd.append("-val")
    if default:
        cmd.append("-default")
    if input:
        cmd.append("-in")
        cmd.append(input)
    if output:
        cmd.append("-out")
        cmd.append(output)
    if solution:
        cmd.append("-sol")
        cmd.append(solution)
    if metric:
        cmd.append("-met")
        cmd.append(metric)
    if A:
        cmd.append("-A")
    if ar:
        cmd.append("-ar")
        cmd.append(str(ar))
    if hausd:
        cmd.append("-hausd")
        cmd.append(str(hausd))
    if hgrad:
        cmd.append("-hgrad")
        cmd.append(str(hgrad))
    if hmax:
        cmd.append("-hmax")
        cmd.append(str(hmax))
    if hmin:
        cmd.append("-hmin")
        cmd.append(str(hmin))
    if hsiz:
        cmd.append("-hsiz")
        cmd.append(str(hsiz))
    if lag or lag == 0:
        if isinstance(lag, bool):
            lag = 0
        cmd.append("-lag")
        cmd.append(str(lag))
    if ls or ls == 0:
        ls_value = ls
        if isinstance(ls_value, bool):
            ls_value = 0
        cmd.append("-ls")
        cmd.append(str(ls_value))
    if _3dMedit:
        cmd.append("-3dMedit")
        cmd.append(str(_3dMedit))

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
        cmd.append("-nreg")
        cmd.append(str(nreg))
    if nsd:
        cmd.append("-nsd")
        cmd.append(str(nsd))
    if optim:
        cmd.append("-optim")
    if opnbdy:
        cmd.append("-opnbdy")
    if rmc:
        cmd.append("-rmc")
        cmd.append(str(rmc))

    subprocess.run(cmd)


def mmgs(
    d=None,
    h=None,
    m=None,
    v=None,
    val=None,
    default=None,
    input=None,
    output=None,
    solution=None,
    metric=None,
    A=None,
    ar=None,
    hausd=None,
    hgrad=None,
    hmax=None,
    hmin=None,
    hsiz=None,
    ls=None,
    noinsert=None,
    nomove=None,
    nosurf=None,
    noswap=None,
    nr=None,
    nreg=None,
    nsd=None,
    optim=None,
    rn=None,
):
    cmd = ["mmgs_O3"]
    if d:
        cmd.append("-d")
    if h:
        cmd.append("-h")
    if m:
        if isinstance(m, bool):
            m = ""
        cmd.append("-m")
        cmd.append(str(m))
    if v:
        if isinstance(v, bool):
            v = 1
        cmd.append("-v")
        cmd.append(str(v))
    if val:
        cmd.append("-val")
    if default:
        cmd.append("-default")
    if input:
        cmd.append("-in")
        cmd.append(input)
    if output:
        cmd.append("-out")
        cmd.append(output)
    if solution:
        cmd.append("-sol")
        cmd.append(solution)
    if metric:
        cmd.append("-met")
        cmd.append(metric)
    if A:
        cmd.append("-A")
    if ar:
        cmd.append("-ar")
        cmd.append(str(ar))
    if hausd:
        cmd.append("-hausd")
        cmd.append(str(hausd))
    if hgrad:
        cmd.append("-hgrad")
        cmd.append(str(hgrad))
    if hmax:
        cmd.append("-hmax")
        cmd.append(str(hmax))
    if hmin:
        cmd.append("-hmin")
        cmd.append(str(hmin))
    if hsiz:
        cmd.append("-hsiz")
        cmd.append(str(hsiz))
    if ls or ls == 0:
        ls_value = ls
        if isinstance(ls_value, bool):
            ls_value = 0
        cmd.append("-ls")
        cmd.append(str(ls_value))
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
        cmd.append("-nreg")
        cmd.append(str(nreg))
    if nsd:
        cmd.append("-nsd")
        cmd.append(str(nsd))
    if optim:
        cmd.append("-optim")
    if rn:
        cmd.append("-rn")

    subprocess.run(cmd)


def mmg3d(
    d=None,
    h=None,
    m=None,
    v=None,
    val=None,
    default=None,
    input=None,
    output=None,
    solution=None,
    metric=None,
    A=None,
    ar=None,
    octree=None,
    hausd=None,
    hgrad=None,
    hmax=None,
    hmin=None,
    hsiz=None,
    lag=None,
    ls=None,
    nofem=None,
    noinsert=None,
    nomove=None,
    nosurf=None,
    noswap=None,
    nr=None,
    nreg=None,
    nsd=None,
    optim=None,
    optimLES=None,
    opnbdy=None,
    rmc=None,
    rn=None,
):
    cmd = ["mmg3d_O3"]
    if d:
        cmd.append("-d")
    if h:
        cmd.append("-h")
    if m:
        if isinstance(m):
            m = ""
        cmd.append("-m")
        cmd.append(str(m))
    if v:
        if isinstance(v):
            v = 1
        cmd.append("-v")
        cmd.append(str(v))
    if val:
        cmd.append("-val")
    if default:
        cmd.append("-default")
    if input:
        cmd.append("-in")
        cmd.append(input)
    if output:
        cmd.append("-out")
        cmd.append(output)
    if solution:
        cmd.append("-sol")
        cmd.append(solution)
    if metric:
        cmd.append("-met")
        cmd.append(metric)
    if A:
        cmd.append("-A")
    if ar:
        cmd.append("-ar")
        cmd.append(str(ar))
    if octree:
        cmd.append("-octree")
        cmd.append(str(octree))
    if hausd:
        cmd.append("-hausd")
        cmd.append(str(hausd))
    if hgrad:
        cmd.append("-hgrad")
        cmd.append(str(hgrad))
    if hmax:
        cmd.append("-hmax")
        cmd.append(str(hmax))
    if hmin:
        cmd.append("-hmin")
        cmd.append(str(hmin))
    if hsiz:
        cmd.append("-hsiz")
        cmd.append(str(hsiz))
    if lag or lag == 0:
        if isinstance(lag, bool):
            lag = 0
        cmd.append("-lag")
        cmd.append(str(lag))
    if ls or ls == 0:
        ls_value = ls
        if isinstance(ls_value, bool):
            ls_value = 0
        cmd.append("-ls")
        cmd.append(str(ls_value))
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
        cmd.append("-nreg")
        cmd.append(str(nreg))
    if nsd:
        cmd.append("-nsd")
        cmd.append(str(nsd))
    if optim:
        cmd.append("-optim")
    if optimLES:
        cmd.append("-optimLES")
    if opnbdy:
        cmd.append("-opnbdy")
    if rmc:
        cmd.append("-rmc")
        cmd.append(str(rmc))
    if rn:
        cmd.append("-rn")

    subprocess.run(cmd)
