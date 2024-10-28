"""
Functions related to external software
    - `Neper - Polycrystal Generation and Meshing`_
    - `mmg - Robust, Open-source & Multidisciplinary Software for Remeshing`_

.. _Neper - Polycrystal Generation and Meshing: https://neper.info/
.. _mmg - Robust, Open-source & Multidisciplinary Software for Remeshing: https://www.mmgtools.org/
"""

from __future__ import annotations

import subprocess

import numpy as np

from microgen.shape import Polyhedron


class Neper:
    @staticmethod
    def run(filename: str, nbCell: int, dimCube: tuple[float, float, float]) -> None:
        """
        Runs neper command from the command line

        command = neper -T -n nbCell -id 1 -dim 3 -domain 'cube(dimCube[0], dimCube[1], dimCube[2])' -morpho gg -o filename

        :param filename: output file
        :param nbCell: Specify the number of cells of the tessellation
        :param dimCube: `neper's documentation`_

        .. _neper's documentation: https://neper.info/doc/neper_t.html#cmdoption-domain
        """

        try:
            # Check if neper is installed and available in the PATH
            subprocess.run(["neper", "--version"], check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print(
                "Neper is not installed. Please install Neper before running this command."
            )
            return

        command = ["neper", "-T", "-n", str(nbCell), "-id", "1", "-dim", "3"]
        command += [
            "-domain",
            "'cube(",
            str(dimCube[0]),
            ",",
            str(dimCube[1]),
            ",",
            str(dimCube[2]),
            ")'",
        ]
        command += ["-morpho", "gg", "-o", filename]
        subprocess.run(command, check=True)

    @staticmethod
    def generateVoronoiFromTessFile(filename: str) -> list[Polyhedron]:
        """
        Generates list of Voronoi polyhedron shapes from a tessellation file generated with neper
        """
        polyhedra = []
        tess = Neper.tessParse(filename)
        for i in range(tess["cells"]["number_of_cells"]):
            global_ind_faces = tess["polyhedra"]["faces"][i]
            global_ind_vertices = []
            for face_ind in global_ind_faces:
                face_vertices = tess["faces"]["vertices"][face_ind - 1]
                global_ind_vertices += face_vertices
            global_ind_vertices = list(
                dict.fromkeys(global_ind_vertices)
            )  # remove duplicates
            global_ind_vertices.sort()

            vertices = []
            local_ind_vertices = []
            for i in range(len(global_ind_vertices)):
                ind = global_ind_vertices[i]
                x = tess["vertices"]["ver_x"][ind - 1]
                y = tess["vertices"]["ver_y"][ind - 1]
                z = tess["vertices"]["ver_z"][ind - 1]
                vertices.append((x, y, z))
                local_ind_vertices.append(i)

            faces = []
            for face_ind in global_ind_faces:
                face_vertices = tess["faces"]["vertices"][face_ind - 1]
                faces.append({"vertices": []})
                for ver in face_vertices:
                    faces[-1]["vertices"].append(global_ind_vertices.index(ver))

            polyhedra.append(Polyhedron(dic={"vertices": vertices, "faces": faces}))
        return polyhedra

    @staticmethod
    def tessParse(filename: str) -> dict[str, dict]:
        """
        Parses tessellation file (.tess) generated with neper.
        Following .tess structure from `neper's`_ documentation:
        Returns a dictionary containing information cells,
        vertices, edges, faces, polyhedra

        .. _neper's: https://neper.info/doc/fileformat.html#tessellation-file-tess
        """
        if len(filename.split(".")) == 1:
            filename += ".tess"
        with open(filename) as f:
            f.readline()  # ***tess
            f.readline()  # **format
            f.readline()  # <format>
            f.readline()  # **general
            f.readline()  # <dim> <type>
            f.readline()  # **cell
            cells = Neper._readCells(f)
            vertices = Neper._readVertices(f)
            edges = Neper._readEdges(f)
            faces = Neper._readFaces(f)
            polyhedra = Neper._readPolyhedra(f)
            # not reading the rest of the file
        return {
            "cells": cells,
            "vertices": vertices,
            "edges": edges,
            "faces": faces,
            "polyhedra": polyhedra,
        }

    @staticmethod
    def _readCells(f):
        cells = {"number_of_cells": int(f.readline())}
        while True:
            tag = f.readline().rstrip()
            if tag == "  *id":
                cells["id"] = []
                while len(cells["id"]) < cells["number_of_cells"]:
                    data = f.readline().split()
                    for id in data:
                        cells["id"].append(int(id))
            elif tag == "  *mode" or tag == "  *modeid":
                cells["mode"] = []
                while len(cells["mode"]) < cells["number_of_cells"]:
                    data = f.readline().split()
                    for id in data:
                        cells["mode"].append(int(id))
            elif tag == "  *seed":
                cells["seed"] = []
                for i in range(cells["number_of_cells"]):
                    data = f.readline().split()
                    cells["seed"].append(
                        {
                            "seed_id": int(data[0]),
                            "seed_x": float(data[1]),
                            "seed_y": float(data[2]),
                            "seed_z": float(data[3]),
                            "seed_weight": float(data[4]),
                        }
                    )
            elif tag == "  *ori":
                cells["ori"] = {"descriptor": f.readline()}
                for i in range(cells["number_of_cells"]):
                    data = f.readline().split()
                    cells["ori"]["cellid_param"] = [
                        float(data[0]),
                        float(data[1]),
                        float(data[2]),
                    ]
            elif tag == "  *orispread":
                pass
            elif tag == "  *lam":
                pass
            elif tag == "  *mode":
                pass
            elif tag == "  *crysym":
                cells["crysym"] = f.readline()
            elif tag == " **vertex":
                return cells
            else:
                raise ValueError(f"tag {tag} not known")

    @staticmethod
    def _readVertices(f):
        vertices = {
            "total_number_of_vertices": int(f.readline()),
            "ver_id": [],
            "ver_x": [],
            "ver_y": [],
            "ver_z": [],
            "ver_state": [],
        }
        for i in range(vertices["total_number_of_vertices"]):
            data = f.readline().split()
            vertices["ver_id"].append(int(data[0]))
            vertices["ver_x"].append(float(data[1]))
            vertices["ver_y"].append(float(data[2]))
            vertices["ver_z"].append(float(data[3]))
            vertices["ver_state"].append(int(data[4]))
        return vertices

    @staticmethod
    def _readEdges(f):
        f.readline()  # **edge
        edges = {
            "total_number_of_edges": int(f.readline().rstrip()),
            "edge_id": [],
            "ver_1": [],
            "ver_2": [],
            "edge_state": [],
        }
        for i in range(edges["total_number_of_edges"]):
            data = f.readline().split()
            edges["edge_id"].append(int(data[0]))
            edges["ver_1"].append(int(data[1]))
            edges["ver_2"].append(int(data[2]))
            edges["edge_state"].append(int(data[3]))
        return edges

    @staticmethod
    def _readFaces(f):
        f.readline()  # **face
        faces = {
            "total_number_of_faces": int(f.readline().rstrip()),
            "face_id": [],
            "vertices": [],
            "edges": [],
        }
        for i in range(faces["total_number_of_faces"]):
            data = f.readline().split()
            faces["face_id"].append(int(data[0]))
            n_ver = int(data[1])
            faces["vertices"].append([])
            for j in range(n_ver):
                faces["vertices"][-1].append(int(data[j + 2]))
            data = f.readline().split()
            n_edg = int(data[0])
            faces["edges"].append([])
            for j in range(n_edg):
                faces["edges"][-1].append(np.abs(int(data[j + 1])))
            f.readline()  # not reading this line
            f.readline()  # not reading this line
        return faces

    @staticmethod
    def _readPolyhedra(f):
        f.readline()  # **polyhedron
        polyhedra = {
            "total_number_of_polyhedra": int(f.readline().rstrip()),
            "poly_id": [],
            "faces": [],
        }
        for i in range(polyhedra["total_number_of_polyhedra"]):
            data = f.readline().split()
            polyhedra["poly_id"] = int(data[0])
            n_fac = int(data[1])
            polyhedra["faces"].append([])
            for j in range(n_fac):
                polyhedra["faces"][-1].append(np.abs(int(data[j + 2])))
        return polyhedra


def parseNeper(filename: str) -> tuple:
    """
    Parses .tess tessellation file obtained with neper
    :param filename: .tess file name
    """

    # First step: read coordinates of seeds, coordinates of vertices, global labels
    # of vertices for each segment, global labels of segments for each face, global
    # labels of faces for each polyhedron from the tessellation file filename.tess
    # generated by NEPER

    file = filename + ".tess"
    fid = open(file)

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

    # Second step: create an object for the polycrystal having a structure similar
    # to the output of compute_voronoi from the py_voro library, compatible with the
    # definition of the polyhedra of the Microgen library

    A = []
    nbPolys = len(polys)
    listeSommets = []
    listeSommetsOut = []
    for i in range(nbPolys):
        voro: dict[str, dict | list] = {}
        voro["original"] = seed[i]
        voro["faces"] = []
        listeSommetsPoly = []
        sommets = []
        for facePoly in polys[i][1:]:
            listeSommetsFace = []
            dicVertices: dict[str, list] = {}
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


class MmgError(Exception):
    """Raised when Mmg command fails"""

    ...


class Mmg:
    @staticmethod
    def _run_mmg_command(cmd: list[str]):
        try:
            subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        except (subprocess.CalledProcessError, FileNotFoundError) as error:
            mmg_failed_command = " ".join(cmd)
            raise MmgError(f"mmg command '{mmg_failed_command}' failed") from error

    @staticmethod
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
            cmd.extend(["-m", str(m)])
        if v:
            if isinstance(v, bool):
                v = 1
            cmd.extend(["-v", str(v)])
        if val:
            cmd.append("-val")
        if default:
            cmd.append("-default")
        if input:
            cmd.extend(["-in", input])
        if output:
            cmd.extend(["-out", output])
        if solution:
            cmd.extend(["-sol", solution])
        if metric:
            cmd.extend(["-met", metric])
        if A:
            cmd.append("-A")
        if ar:
            cmd.extend(["-ar", str(ar)])
        if hausd:
            cmd.extend(["-hausd", str(hausd)])
        if hgrad:
            cmd.extend(["-hgrad", str(hgrad)])
        if hmax:
            cmd.extend(["-hmax", str(hmax)])
        if hmin:
            cmd.extend(["-hmin", str(hmin)])
        if hsiz:
            cmd.extend(["-hsiz", str(hsiz)])
        if lag or lag == 0:
            if isinstance(lag, bool):
                lag = 0
            cmd.extend(["-lag", str(lag)])
        if ls or ls == 0:
            ls_value = ls
            if isinstance(ls_value, bool):
                ls_value = 0
            cmd.extend(["-ls", str(ls_value)])
        if _3dMedit:
            cmd.extend(["-3dMedit", str(_3dMedit)])

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
            cmd.extend(["-nreg", str(nreg)])
        if nsd:
            cmd.extend(["-nsd", str(nsd)])
        if optim:
            cmd.append("-optim")
        if opnbdy:
            cmd.append("-opnbdy")
        if rmc:
            cmd.extend(["-rmc", str(rmc)])

        Mmg._run_mmg_command(cmd)

    @staticmethod
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
            cmd.extend(["-m", str(m)])
        if v:
            if isinstance(v, bool):
                v = 1
            cmd.extend(["-v", str(v)])
        if val:
            cmd.append("-val")
        if default:
            cmd.append("-default")
        if input:
            cmd.extend(["-in", input])
        if output:
            cmd.extend(["-out", output])
        if solution:
            cmd.extend(["-sol", solution])
        if metric:
            cmd.extend(["-met", metric])
        if A:
            cmd.append("-A")
        if ar:
            cmd.extend(["-ar", str(ar)])
        if hausd:
            cmd.extend(["-hausd", str(hausd)])
        if hgrad:
            cmd.extend(["-hgrad", str(hgrad)])
        if hmax:
            cmd.extend(["-hmax", str(hmax)])
        if hmin:
            cmd.extend(["-hmin", str(hmin)])
        if hsiz:
            cmd.extend(["-hsiz", str(hsiz)])
        if ls or ls == 0:
            ls_value = ls
            if isinstance(ls_value, bool):
                ls_value = 0
            cmd.extend(["-ls", str(ls_value)])
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
            cmd.extend(["-nreg", str(nreg)])
        if nsd:
            cmd.extend(["-nsd", str(nsd)])
        if optim:
            cmd.append("-optim")
        if rn:
            cmd.append("-rn")
        Mmg._run_mmg_command(cmd)

    @staticmethod
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
            if isinstance(m, bool):
                m = ""
            cmd.extend(["-m", str(m)])
        if v:
            if isinstance(v, bool):
                v = 1
            cmd.extend(["-v", str(v)])
        if val:
            cmd.append("-val")
        if default:
            cmd.append("-default")
        if input:
            cmd.extend(["-in", input])
        if output:
            cmd.extend(["-out", output])
        if solution:
            cmd.extend(["-sol", solution])
        if metric:
            cmd.extend(["-met", metric])
        if A:
            cmd.append("-A")
        if ar:
            cmd.extend(["-ar", str(ar)])
        if octree:
            cmd.extend(["-octree", str(octree)])
        if hausd:
            cmd.extend(["-hausd", str(hausd)])
        if hgrad:
            cmd.extend(["-hgrad", str(hgrad)])
        if hmax:
            cmd.extend(["-hmax", str(hmax)])
        if hmin:
            cmd.extend(["-hmin", str(hmin)])
        if hsiz:
            cmd.extend(["-hsiz", str(hsiz)])
        if lag or lag == 0:
            if isinstance(lag, bool):
                lag = 0
            cmd.extend(["-lag", str(lag)])
        if ls or ls == 0:
            ls_value = ls
            if isinstance(ls_value, bool):
                ls_value = 0
            cmd.extend(["-ls", str(ls_value)])
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
            cmd.extend(["-nreg", str(nreg)])
        if nsd:
            cmd.extend(["-nsd", str(nsd)])
        if optim:
            cmd.append("-optim")
        if optimLES:
            cmd.append("-optimLES")
        if opnbdy:
            cmd.append("-opnbdy")
        if rmc:
            cmd.extend(["-rmc", str(rmc)])
        if rn:
            cmd.append("-rn")

        Mmg._run_mmg_command(cmd)
