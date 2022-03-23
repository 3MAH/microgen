import numpy as np
import gmsh

gmsh.initialize()


def mesh(mesh_file, listPhases, size, order=1):

    flatListPhases = [phase for phase_list in listPhases for phase in phase_list]
    nbTags = len(flatListPhases)
    listTagsNb = [len(phase_list) for phase_list in listPhases]

    print(nbTags)
    print(listTagsNb)
    FlatListTags = list(range(1, nbTags + 1, 1))
    print(FlatListTags)

    listTags = []
    index = 0
    for i, phase_list in enumerate(listPhases):
        temp = []
        for j, phase in enumerate(phase_list):
            index = index + 1
            temp.append(index)
        listTags.append(temp)
    print(listTags)

    listDimTags = [(3, tag) for tag in FlatListTags]
    print(listDimTags)

    ToMesh = gmsh.model.occ.importShapes(mesh_file, highestDimOnly=True)
    print(ToMesh)

    #    #get all elementary entities in the model
    #    entities = gmsh.model.getEntities()
    #    print(len(entities))
    #
    #    for e in entities:
    #        print("Entity " + str(e) + " of type " + gmsh.model.getType(e[0], e[1]))
    #        # get the mesh nodes for each elementary entity
    #        nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(e[0], e[1])
    #        # get the mesh elements for each elementary entity
    #        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[0], e[1])
    #        # count number of elements
    #        numElem = sum(len(i) for i in elemTags)
    #        print(" - mesh has " + str(len(nodeTags)) + " nodes and " + str(numElem) +
    #              " elements")
    #        boundary = gmsh.model.getBoundary([e])
    #        print(" - boundary entities " + str(boundary))
    #        partitions = gmsh.model.getPartitions(e[0], e[1])
    #        if len(partitions):
    #            print(" - Partition tag(s): " + str(partitions) + " - parent entity " +
    #                  str(gmsh.model.getParent(e[0], e[1])))
    #        for t in elemTypes:
    #            name, dim, order, numv, parv, _ = gmsh.model.mesh.getElementProperties(
    #                t)
    #            print(" - Element type: " + name + ", order " + str(order) + " (" +
    #                  str(numv) + " nodes in param coord: " + str(parv) + ")")

    if len(listDimTags) > 1:
        print(listDimTags[:-1])
        print(listDimTags[-1])
        outDimTags, outDimTagsMap = gmsh.model.occ.fragment(
            listDimTags[:-1], [listDimTags[-1]]
        )

    #    ent = gmsh.model.getEntities(3)
    #    print(ent)
    #    gmsh.model.occ.healShapes()
    gmsh.model.occ.synchronize()

    #    print(a)
    #    print(len(a))
    #    gmsh.model.occ.fuse([(3, 3)], [(3, 4)])
    #    gmsh.model.occ.fuse([(3, 5)], [(3, 6)])

    for i, tag in enumerate(listTags):
        print(i)
        print(type(i))
        ps_i = gmsh.model.addPhysicalGroup(3, tag)
        gmsh.model.setPhysicalName(3, ps_i, "Mat" + str(i))

    # gmsh.model.mesh.setSize(gmsh.model.getEntities(0), Size)
    # p = gmsh.model.getBoundary(ToMesh, False, False, True)  # Get all points
    p = gmsh.model.getEntities()
    #    print(p)

    gmsh.model.mesh.setSize(p, size)
    gmsh.model.mesh.setOrder(order)
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.write("Mesh.msh")
    gmsh.finalize()


def meshPeriodic(mesh_file, rve, listPhases, size, order=1):

    flatListPhases = [phase for phase_list in listPhases for phase in phase_list]
    nbTags = len(flatListPhases)
    listTagsNb = [len(phase_list) for phase_list in listPhases]

    print(nbTags)
    print(listTagsNb)
    flatListTags = list(range(1, nbTags + 1, 1))
    print(flatListTags)

    listTags = []
    index = 0
    for i, phase_list in enumerate(listPhases):
        temp = []
        for j, phase in enumerate(phase_list):
            index = index + 1
            temp.append(index)
        listTags.append(temp)
    print(listTags)

    listDimTags = [(3, tag) for tag in flatListTags]
    print(listDimTags)

    toMesh = gmsh.model.occ.importShapes(mesh_file, highestDimOnly=True)
    print(toMesh)

    #    #get all elementary entities in the model
    #    entities = gmsh.model.getEntities()
    #    print(len(entities))
    #
    #    for e in entities:
    #        print("Entity " + str(e) + " of type " + gmsh.model.getType(e[0], e[1]))
    #        # get the mesh nodes for each elementary entity
    #        nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(e[0], e[1])
    #        # get the mesh elements for each elementary entity
    #        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[0], e[1])
    #        # count number of elements
    #        numElem = sum(len(i) for i in elemTags)
    #        print(" - mesh has " + str(len(nodeTags)) + " nodes and " + str(numElem) +
    #              " elements")
    #        boundary = gmsh.model.getBoundary([e])
    #        print(" - boundary entities " + str(boundary))
    #        partitions = gmsh.model.getPartitions(e[0], e[1])
    #        if len(partitions):
    #            print(" - Partition tag(s): " + str(partitions) + " - parent entity " +
    #                  str(gmsh.model.getParent(e[0], e[1])))
    #        for t in elemTypes:
    #            name, dim, order, numv, parv, _ = gmsh.model.mesh.getElementProperties(
    #                t)
    #            print(" - Element type: " + name + ", order " + str(order) + " (" +
    #                  str(numv) + " nodes in param coord: " + str(parv) + ")")

    size_box = np.min(np.array([rve.dx, rve.dy, rve.dz]))
    eps = 1.0e-3 * size_box
    if len(listDimTags) > 1:
        outDimTags, outDimTagsMap = gmsh.model.occ.fragment(
            listDimTags[:-1], [listDimTags[-1]]
        )
    gmsh.model.occ.synchronize()

    for i, tag in enumerate(listTags):
        print(i)
        print(type(i))
        ps_i = gmsh.model.addPhysicalGroup(3, tag)
        gmsh.model.setPhysicalName(3, ps_i, "Mat" + str(i))

    #    p = gmsh.model.getBoundary(ToMesh, False, False, True)  # Get all points
    p = gmsh.model.getEntities()
    #    print(p)

    # We now identify corresponding surfaces on the left and right sides of the
    # geometry automatically.

    # We get all the entities on the Xm
    translation = [1, 0, 0, rve.dx, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    sxmin = gmsh.model.getEntitiesInBoundingBox(
        0 - eps, -eps, -eps, eps, rve.dy + eps, rve.dy + eps, 2
    )

    for i in sxmin:
        # Then we get the bounding box of each left surface
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(i[0], i[1])
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
        for j in sxmax:
            xmin2, ymin2, zmin2, xmax2, ymax2, zmax2 = gmsh.model.getBoundingBox(
                j[0], j[1]
            )
            xmin2 -= 1
            xmax2 -= 1

            print(xmin, ymin, zmin, xmax, ymax, zmax)
            # ...and if they match, we apply the periodicity constraint
            if (
                abs(xmin2 - xmin) < eps
                and abs(xmax2 - xmax) < eps
                and abs(ymin2 - ymin) < eps
                and abs(ymax2 - ymax) < eps
                and abs(zmin2 - zmin) < eps
                and abs(zmax2 - zmax) < eps
            ):
                gmsh.model.mesh.setPeriodic(2, [j[1]], [i[1]], translation)

    # We get all the entities on the Ym
    translation = [1, 0, 0, 0, 0, 1, 0, rve.dy, 0, 0, 1, 0, 0, 0, 0, 1]
    symin = gmsh.model.getEntitiesInBoundingBox(
        0 - eps, -eps, -eps, rve.dx + eps, eps, rve.dz + eps, 2
    )

    for i in symin:
        # Then we get the bounding box of each left surface
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(i[0], i[1])
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
        for j in symax:
            xmin2, ymin2, zmin2, xmax2, ymax2, zmax2 = gmsh.model.getBoundingBox(
                j[0], j[1]
            )
            ymin2 -= 1
            ymax2 -= 1

            print(xmin, ymin, zmin, xmax, ymax, zmax)
            # ...and if they match, we apply the periodicity constraint
            if (
                abs(xmin2 - xmin) < eps
                and abs(xmax2 - xmax) < eps
                and abs(ymin2 - ymin) < eps
                and abs(ymax2 - ymax) < eps
                and abs(zmin2 - zmin) < eps
                and abs(zmax2 - zmax) < eps
            ):
                gmsh.model.mesh.setPeriodic(2, [j[1]], [i[1]], translation)

    # We get all the entities on the Zm
    translation = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, rve.dz, 0, 0, 0, 1]
    szmin = gmsh.model.getEntitiesInBoundingBox(
        0 - eps, -eps, -eps, rve.dx + eps, rve.dy + eps, eps, 2
    )

    for i in szmin:
        # Then we get the bounding box of each left surface
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(i[0], i[1])
        # We translate the bounding box to the right and look for surfaces inside
        # it:
        print(xmin, ymin, zmin, xmax, ymax, zmax)
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
        for j in szmax:
            xmin2, ymin2, zmin2, xmax2, ymax2, zmax2 = gmsh.model.getBoundingBox(
                j[0], j[1]
            )
            zmin2 -= 1
            zmax2 -= 1

            print(xmin, ymin, zmin, xmax, ymax, zmax)
            # ...and if they match, we apply the periodicity constraint
            if (
                abs(xmin2 - xmin) < eps
                and abs(xmax2 - xmax) < eps
                and abs(ymin2 - ymin) < eps
                and abs(ymax2 - ymax) < eps
                and abs(zmin2 - zmin) < eps
                and abs(zmax2 - zmax) < eps
            ):
                gmsh.model.mesh.setPeriodic(2, [j[1]], [i[1]], translation)

    p = gmsh.model.getEntities()
    gmsh.model.mesh.setSize(p, size)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.write("MeshPeriodic.msh")
    gmsh.finalize()
