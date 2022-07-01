"""
Mesh using gmsh
"""

import gmsh
import numpy as np

from .phase import Phase
from .rve import Rve


def mesh(
    mesh_file: str,
    listPhases: list[Phase],
    size: float,
    order: int,
    output_file: str = "Mesh.msh",
    mshFileVersion: int = 4,
) -> None:
    """
    Meshes step file with gmsh with list of phases management

    :param mesh_file: step file to mesh
    :param listPhases: list of phases to mesh
    :param size: mesh size constraint (see: `gmsh.model.mesh.setSize(dimTags, size)`_)
    :param order: see `gmsh.model.mesh.setOrder(order)`_
    :param output_file: output file (.msh, .vtk)
    :param mshFileVersion: gmsh file version

    .. _gmsh.model.mesh.setOrder(order): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L1688
    .. _gmsh.model.mesh.setSize(dimTags, size): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L3140
    """
    gmsh.initialize()
    gmsh.option.setNumber(
        name="General.Verbosity", value=1
    )  # this would still print errors, but not warnings

    gmsh.model.mesh.setOrder(order=order)
    gmsh.option.setNumber(name="Mesh.MshFileVersion", value=mshFileVersion)

    flatListSolids = [solid for phase in listPhases for solid in phase.solids]
    nbTags = len(flatListSolids)
    FlatListTags = list(range(1, nbTags + 1, 1))

    listTags = []
    index = 0
    for i, phase in enumerate(listPhases):
        temp = []
        for j, solid in enumerate(phase.solids):
            index = index + 1
            temp.append(index)
        listTags.append(temp)

    listDimTags = [(3, tag) for tag in FlatListTags]

    gmsh.model.occ.importShapes(fileName=mesh_file, highestDimOnly=True)

    if len(listDimTags) > 1:
        gmsh.model.occ.fragment(
            objectDimTags=listDimTags[:-1], toolDimTags=[listDimTags[-1]]
        )

    gmsh.model.occ.synchronize()

    for i, tag in enumerate(listTags):
        ps_i = gmsh.model.addPhysicalGroup(dim=3, tags=tag)
        gmsh.model.setPhysicalName(dim=3, tag=ps_i, name="Mat" + str(i))

    p = gmsh.model.getEntities()

    gmsh.model.mesh.setSize(dimTags=p, size=size)
    gmsh.model.mesh.generate(dim=3)
    gmsh.write(fileName=output_file)
    gmsh.finalize()


def meshPeriodic(
    mesh_file: str,
    rve: Rve,
    listPhases: list[Phase],
    size: float,
    order: int,
    output_file: str = "MeshPeriodic.msh",
    mshFileVersion: int = 4,
) -> None:
    """
    Meshes periodic geometries with gmsh

    :param mesh_file: step file to mesh
    :param rve: RVE for periodicity
    :param listPhases: list of phases to mesh
    :param size: mesh size constraint (see: `gmsh.model.mesh.setSize(dimTags, size)`_)
    :param order: see `gmsh.model.mesh.setOrder(order)`_
    :param output_file: output file (.msh, .vtk)
    :param mshFileVersion: gmsh file version

    .. _gmsh.model.mesh.setOrder(order): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L1688
    .. _gmsh.model.mesh.setSize(dimTags, size): https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py#L3140
    """
    gmsh.initialize()
    gmsh.option.setNumber(
        "General.Verbosity", 1
    )  # this would still print errors, but not warnings

    gmsh.model.mesh.setOrder(order)
    gmsh.option.setNumber("Mesh.MshFileVersion", mshFileVersion)

    flatListSolids = [solid for phase in listPhases for solid in phase.solids]
    nbTags = len(flatListSolids)

    flatListTags = list(range(1, nbTags + 1, 1))

    listTags = []
    index = 0
    for i, phase in enumerate(listPhases):
        temp = []
        for j, solid in enumerate(phase.solids):
            index = index + 1
            temp.append(index)
        listTags.append(temp)

    listDimTags = [(3, tag) for tag in flatListTags]

    gmsh.model.occ.importShapes(mesh_file, highestDimOnly=True)

    size_box = np.min(np.array([rve.dx, rve.dy, rve.dz]))
    eps = 1.0e-3 * size_box
    if len(listDimTags) > 1:
        outDimTags, outDimTagsMap = gmsh.model.occ.fragment(
            listDimTags[:-1], [listDimTags[-1]]
        )
    gmsh.model.occ.synchronize()

    for i, tag in enumerate(listTags):
        ps_i = gmsh.model.addPhysicalGroup(3, tag)
        gmsh.model.setPhysicalName(3, ps_i, "Mat" + str(i))

    p = gmsh.model.getEntities()

    # We now identify corresponding surfaces on the left and right sides of the
    # geometry automatically.

    # We get all the entities on the Xm
    translation = [1, 0, 0, rve.dx, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    sxmin = gmsh.model.getEntitiesInBoundingBox(
        0 - eps, -eps, -eps, eps, rve.dy + eps, rve.dy + eps, 2
    )

    for tup_min in sxmin:
        # Then we get the bounding box of each left surface
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(
            tup_min[0], tup_min[1]
        )
        # We translate the bounding box to the right and look for surfaces inside
        # it:
        sxmax = gmsh.model.getEntitiesInBoundingBox(
            xmin - eps + 1,
            ymin - eps,
            zmin - eps,
            xmax + eps + 1,
            ymax + eps,
            zmax + eps,
            2,
        )
        # For all the matches, we compare the corresponding bounding boxes...
        for tup_max in sxmax:
            xmin2, ymin2, zmin2, xmax2, ymax2, zmax2 = gmsh.model.getBoundingBox(
                tup_max[0], tup_max[1]
            )
            xmin2 -= 1
            xmax2 -= 1

            # ...and if they match, we apply the periodicity constraint
            if (
                abs(xmin2 - xmin) < eps
                and abs(xmax2 - xmax) < eps
                and abs(ymin2 - ymin) < eps
                and abs(ymax2 - ymax) < eps
                and abs(zmin2 - zmin) < eps
                and abs(zmax2 - zmax) < eps
            ):
                gmsh.model.mesh.setPeriodic(2, [tup_max[1]], [tup_min[1]], translation)

    # We get all the entities on the Ym
    translation = [1, 0, 0, 0, 0, 1, 0, rve.dy, 0, 0, 1, 0, 0, 0, 0, 1]
    symin = gmsh.model.getEntitiesInBoundingBox(
        0 - eps, -eps, -eps, rve.dx + eps, eps, rve.dz + eps, 2
    )

    for tup_min in symin:
        # Then we get the bounding box of each left surface
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(
            tup_min[0], tup_min[1]
        )
        # We translate the bounding box to the right and look for surfaces inside
        # it:
        symax = gmsh.model.getEntitiesInBoundingBox(
            xmin - eps,
            ymin - eps + 1,
            zmin - eps,
            xmax + eps,
            ymax + eps + 1,
            zmax + eps,
            2,
        )
        # For all the matches, we compare the corresponding bounding boxes...
        for tup_max in symax:
            xmin2, ymin2, zmin2, xmax2, ymax2, zmax2 = gmsh.model.getBoundingBox(
                tup_max[0], tup_max[1]
            )
            ymin2 -= 1
            ymax2 -= 1

            # ...and if they match, we apply the periodicity constraint
            if (
                abs(xmin2 - xmin) < eps
                and abs(xmax2 - xmax) < eps
                and abs(ymin2 - ymin) < eps
                and abs(ymax2 - ymax) < eps
                and abs(zmin2 - zmin) < eps
                and abs(zmax2 - zmax) < eps
            ):
                gmsh.model.mesh.setPeriodic(2, [tup_max[1]], [tup_min[1]], translation)

    # We get all the entities on the Zm
    translation = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, rve.dz, 0, 0, 0, 1]
    szmin = gmsh.model.getEntitiesInBoundingBox(
        0 - eps, -eps, -eps, rve.dx + eps, rve.dy + eps, eps, 2
    )

    for tup_min in szmin:
        # Then we get the bounding box of each left surface
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(
            tup_min[0], tup_min[1]
        )
        # We translate the bounding box to the right and look for surfaces inside
        # it:
        szmax = gmsh.model.getEntitiesInBoundingBox(
            xmin - eps,
            ymin - eps,
            zmin - eps + 1,
            xmax + eps,
            ymax + eps,
            zmax + eps + 1,
            2,
        )
        # For all the matches, we compare the corresponding bounding boxes...
        for tup_max in szmax:
            xmin2, ymin2, zmin2, xmax2, ymax2, zmax2 = gmsh.model.getBoundingBox(
                tup_max[0], tup_max[1]
            )
            zmin2 -= 1
            zmax2 -= 1

            # ...and if they match, we apply the periodicity constraint
            if (
                abs(xmin2 - xmin) < eps
                and abs(xmax2 - xmax) < eps
                and abs(ymin2 - ymin) < eps
                and abs(ymax2 - ymax) < eps
                and abs(zmin2 - zmin) < eps
                and abs(zmax2 - zmax) < eps
            ):
                gmsh.model.mesh.setPeriodic(2, [tup_max[1]], [tup_min[1]], translation)

    p = gmsh.model.getEntities()
    gmsh.model.mesh.setSize(p, size)
    gmsh.model.mesh.generate(3)
    gmsh.write(output_file)
    gmsh.finalize()
