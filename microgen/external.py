"""
Functions related to external softwares
    - `Neper - Polycrystal Generation and Meshing`_
    - `mmg - Robust, Open-source & Multidisciplinary Software for Remeshing`_

.. _Neper - Polycrystal Generation and Meshing: https://neper.info/
.. _mmg - Robust, Open-source & Multidisciplinary Software for Remeshing: https://www.mmgtools.org/
"""

import subprocess
from typing import Union

import numpy as np


def lanceNeper(filename: str, nbCell: int, dimCube: list) -> None:
    """
    Runs neper command from the command line
    """
    command = "neper -T -n " + str(nbCell) + " -id 1 -dim 3"
    command = (
        command
        + " -domain 'cube("
        + str(dimCube[0])
        + ","
        + str(dimCube[1])
        + ","
        + str(dimCube[2])
        + ")'"
    )
    command = command + " -morpho " "gg" " -o " + filename

    try:
        subprocess.check_output(command.split(" "), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print(
            "neper command did not work, check if it is installed or contact a developer"
        )


def parseNeper(filename: str) -> tuple:
    """
    Parses file obtained with neper
    """

    # 1ere étape : lecture des coordonées des germes, coordonées des sommets,
    # labels globaux des sommets pour chaque segment, labels globaux des segments
    # pour chaque face, labels globaux des faces pour chaque polyèdre à partir
    # du fichier de tesselation filename.tess généré par NEPER

    file = filename + ".tess"
    fid = open(file, "r")

    flagChercheSeed = False
    seed = []

    flagChercheVertex = False
    vertices = []

    flagChercheEdges = False
    edges = []

    flagChercheFaces = False
    faces = []
    iLigneFace = 1

    flagCherchePolys = False
    polys = []

    flagExtraction = False
    i = 0

    for line in fid:

        if flagExtraction:
            if "*" in line:
                flagExtraction = False
            else:
                if i == 1:  # Extraction des germes
                    ligneDecoupee = line.split()
                    seed.append(
                        [
                            float(ligneDecoupee[1]),
                            float(ligneDecoupee[2]),
                            float(ligneDecoupee[3]),
                        ]
                    )
                if i == 2:  # Extraction des sommets
                    ligneDecoupee = line.split()
                    vertices.append(
                        [
                            int(ligneDecoupee[0]),
                            float(ligneDecoupee[1]),
                            float(ligneDecoupee[2]),
                            float(ligneDecoupee[3]),
                        ]
                    )
                if i == 3:  # Extraction des segments
                    ligneDecoupee = line.split()
                    edges.append(
                        [
                            int(ligneDecoupee[0]),
                            int(ligneDecoupee[1]),
                            int(ligneDecoupee[2]),
                        ]
                    )
                if i == 4:  # Extraction des faces
                    if iLigneFace == 1:
                        tmp = []
                        tmp.append(int(line.split()[0]))
                    if iLigneFace == 2:
                        ligneDecoupee = line.split()
                        nbEdges = int(ligneDecoupee[0])
                        for k in range(nbEdges):
                            tmp.append(int(ligneDecoupee[k + 1]))
                        faces.append(tmp)
                    if iLigneFace == 4:
                        iLigneFace = 0

                    iLigneFace += 1
                if i == 5:  # Extraction des cellules
                    tmp = []
                    ligneDecoupee = line.split()
                    tmp.append(int(ligneDecoupee[0]))
                    nbFaces = int(ligneDecoupee[1])
                    for k in range(nbFaces):
                        tmp.append(np.abs(int(ligneDecoupee[k + 2])))
                    polys.append(tmp)

        if "*seed" in line:
            flagChercheSeed = True
            i += 1

        if "**vertex" in line:
            flagChercheVertex = True
            i += 1
            next(fid)

        if "**edge" in line:
            flagChercheEdges = True
            i += 1
            next(fid)

        if "**face" in line:
            flagChercheFaces = True
            i += 1
            next(fid)

        if "**polyhedron" in line:
            flagCherchePolys = True
            i += 1
            next(fid)

        if flagChercheSeed:
            flagExtraction = True
            flagChercheSeed = False

        if flagChercheVertex:
            flagExtraction = True
            flagChercheVertex = False

        if flagChercheEdges:
            flagExtraction = True
            flagChercheEdges = False

        if flagChercheFaces:
            flagExtraction = True
            flagChercheFaces = False

        if flagCherchePolys:
            flagExtraction = True
            flagCherchePolys = False

    fid.close()

    # 2e étape : création d'un objet pour le polycristal ayant une structure semblable
    # à la sortie de compute_voronoi de la librairie py_voro, compatible avec la définition
    # des polyèdres de la librairie Microgen

    A = []
    nbPolys = len(polys)
    listeSommets = []
    listeSommetsOut = []
    for i in range(nbPolys):
        voro = {}  # type: dict[str, Union[dict, list]]
        voro["original"] = seed[i]
        voro["faces"] = []
        listeSommetsPoly = []
        sommets = []
        for facePoly in polys[i][1:]:
            listeSommetsFace = []
            dicVertices = {}  # type: dict[str, list]
            dicVertices["vertices"] = []
            for segment in faces[facePoly - 1][1:]:
                # Listes des sommets avec numérotation globale et faces associées
                if edges[np.abs(segment) - 1][1] not in listeSommets:
                    dicFacesAssociees = {}
                    dicFacesAssociees["cell_associees"] = [facePoly]
                    listeSommets.append(edges[np.abs(segment) - 1][1])
                    listeSommetsOut.append(
                        (edges[np.abs(segment) - 1][1], dicFacesAssociees)
                    )
                else:
                    idx = listeSommets.index(edges[np.abs(segment) - 1][1])
                    if facePoly not in listeSommetsOut[idx][1]["cell_associees"]:
                        listeSommetsOut[idx][1]["cell_associees"].append(facePoly)
                    # sommets.append(vertices[edges[np.abs(segment)-1][1]-1][1:])
                if edges[np.abs(segment) - 1][2] not in listeSommets:
                    dicFacesAssociees = {}
                    dicFacesAssociees["cell_associees"] = [facePoly]
                    listeSommets.append(edges[np.abs(segment) - 1][2])
                    listeSommetsOut.append(
                        (edges[np.abs(segment) - 1][2], dicFacesAssociees)
                    )
                else:
                    idx = listeSommets.index(edges[np.abs(segment) - 1][2])
                    if facePoly not in listeSommetsOut[idx][1]["cell_associees"]:
                        listeSommetsOut[idx][1]["cell_associees"].append(facePoly)
                    # sommets.append(vertices[edges[np.abs(segment)-1][2]-1][1:])

                # Listes des sommets et de leurs coordonnées avec numérotation locale
                # à inclure dans voro
                if segment < 0:
                    if edges[np.abs(segment) - 1][2] not in listeSommetsPoly:
                        listeSommetsPoly.append(edges[np.abs(segment) - 1][2])
                        sommets.append(vertices[edges[np.abs(segment) - 1][2] - 1][1:])
                    if edges[np.abs(segment) - 1][2] not in listeSommetsFace:
                        listeSommetsFace.append(edges[np.abs(segment) - 1][2])
                        dicVertices["vertices"].append(
                            listeSommetsPoly.index(edges[np.abs(segment) - 1][2])
                        )

                    if edges[np.abs(segment) - 1][1] not in listeSommetsPoly:
                        listeSommetsPoly.append(edges[np.abs(segment) - 1][1])
                        sommets.append(vertices[edges[np.abs(segment) - 1][1] - 1][1:])
                    if edges[np.abs(segment) - 1][1] not in listeSommetsFace:
                        listeSommetsFace.append(edges[np.abs(segment) - 1][1])
                        dicVertices["vertices"].append(
                            listeSommetsPoly.index(edges[np.abs(segment) - 1][1])
                        )

                if segment > 0:
                    if edges[np.abs(segment) - 1][1] not in listeSommetsPoly:
                        listeSommetsPoly.append(edges[np.abs(segment) - 1][1])
                        sommets.append(vertices[edges[np.abs(segment) - 1][1] - 1][1:])
                    if edges[np.abs(segment) - 1][1] not in listeSommetsFace:
                        listeSommetsFace.append(edges[np.abs(segment) - 1][1])
                        dicVertices["vertices"].append(
                            listeSommetsPoly.index(edges[np.abs(segment) - 1][1])
                        )

                    if edges[np.abs(segment) - 1][2] not in listeSommetsPoly:
                        listeSommetsPoly.append(edges[np.abs(segment) - 1][2])
                        sommets.append(vertices[edges[np.abs(segment) - 1][2] - 1][1:])
                    if edges[np.abs(segment) - 1][2] not in listeSommetsFace:
                        listeSommetsFace.append(edges[np.abs(segment) - 1][2])
                        dicVertices["vertices"].append(
                            listeSommetsPoly.index(edges[np.abs(segment) - 1][2])
                        )
            voro["faces"].append(dicVertices)

        voro["vertices"] = sommets
        A.append(voro)

    return A, seed, listeSommetsOut, edges, faces, polys


# fin ajout


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
    """
    Runs mmg2d_O3 from the command line with given arguments
    """
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

    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print(
            "mmg command did not work, check if it is installed or contact a developer"
        )


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
    """
    Runs mmgs_O3 from the command line with given arguments
    """
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

    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print(
            "mmg command did not work, check if it is installed or contact a developer"
        )


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
    """
    Runs mmg3d_O3 from the command line with given arguments
    """
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

    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print(
            "mmg command did not work, check if it is installed or contact a developer"
        )
