import gmsh
import numpy as np
import math as m

# Remeshing while preserving mesh periodicity

def remesh_periodic(mesh_file: str, rve: Rve, output_file: str = "mesh_reqtri.mesh"
) -> None:
    """
    Prepares the mesh file for periodic remeshing

    :param mesh_file: mesh file to remesh
    :param rve: RVE for periodicity
    :param output_file: output file (must be .mesh)
    """

    # Read mesh

    gmsh.initialize()
    gmsh.open(mesh_file)

    ntets = len(gmsh.model.mesh.getElements(3)[1][0])

    gmsh.model.mesh.createTopology() # Creates boundary entities

    boundarySurface = gmsh.model.getEntities(2)
    boundaryElementsType, boundaryElementsTags, boundaryElementsNodes = gmsh.model.mesh.getElements(2)
    boundaryElementsNodesCoords = np.zeros((len(boundaryElementsNodes[0]), 3))

    for i in range(len(boundaryElementsNodes[0])):
        boundaryElementsNodesCoords[i] = gmsh.model.mesh.getNode(boundaryElementsNodes[0][i])[0]

    # Find and keep only boundary elements (here, triangles) that are on the boundary of the rve

    def computeNormal(node1Coords, node2Coords, node3Coords):
        """Computes normal vector of a triangle in mesh, defined by the coordinates of its three vertices"""
        u = node2Coords-node1Coords
        v = node3Coords-node1Coords
        n = np.cross(u,v)
        normal = n/np.linalg.norm(n)
        return normal

    def isxmin(node1Coords, node2Coords, node3Coords, rve):
        """Determines whether a triangle (defined by its 3 nodes) is on xmin boundary of a rve"""
        xmin = rve.xmin
        xminNormal = np.array([-1.0, 0.0, 0.0])
        if m.isclose(node1Coords[0], xmin) and m.isclose(node2Coords[0], xmin) and m.isclose(node3Coords[0], xmin) and np.allclose(computeNormal(node1Coords, node2Coords, node3Coords), xminNormal):
            return True
        return False

    def isxmax(node1Coords, node2Coords, node3Coords, rve):
        """Determines whether a triangle (defined by its 3 nodes) is on xmax boundary of a rve"""
        xmax = rve.xmax
        xmaxNormal = np.array([1.0, 0.0, 0.0])
        if m.isclose(node1Coords[0], xmax) and m.isclose(node2Coords[0], xmax) and m.isclose(node3Coords[0], xmax) and np.allclose(computeNormal(node1Coords, node2Coords, node3Coords), xmaxNormal):
            return True
        return False

    def isymin(node1Coords, node2Coords, node3Coords, rve):
        """Determines whether a triangle (defined by its 3 nodes) is on ymin boundary of a rve"""
        ymin = rve.ymin
        yminNormal = np.array([0.0, -1.0, 0.0])
        if m.isclose(node1Coords[1], ymin) and m.isclose(node2Coords[1], ymin) and m.isclose(node3Coords[1], ymin) and np.allclose(computeNormal(node1Coords, node2Coords, node3Coords), yminNormal):
            return True
        return False

    def isymax(node1Coords, node2Coords, node3Coords, rve):
        """Determines whether a triangle (defined by its 3 nodes) is on ymax boundary of a rve"""
        ymax = rve.ymax
        ymaxNormal = np.array([0.0, 1.0, 0.0])
        if m.isclose(node1Coords[1], ymax) and m.isclose(node2Coords[1], ymax) and m.isclose(node3Coords[1], ymax) and np.allclose(computeNormal(node1Coords, node2Coords, node3Coords), ymaxNormal):
            return True
        return False

    def iszmin(node1Coords, node2Coords, node3Coords, rve):
        """Determines whether a triangle (defined by its 3 nodes) is on zmin boundary of a rve"""
        zmin = rve.zmin
        zminNormal = np.array([0.0, 0.0, -1.0])
        if m.isclose(node1Coords[2], zmin) and m.isclose(node2Coords[2], zmin) and m.isclose(node3Coords[2], zmin) and np.allclose(computeNormal(node1Coords, node2Coords, node3Coords), zminNormal):
            return True
        return False

    def iszmax(node1Coords, node2Coords, node3Coords, rve):
        """Determines whether a triangle (defined by its 3 nodes) is on zmax boundary of a rve"""
        zmax = rve.zmax
        zmaxNormal = np.array([0.0, 0.0, 1.0])
        if m.isclose(node1Coords[2], zmax) and m.isclose(node2Coords[2], zmax) and m.isclose(node3Coords[2], zmax) and np.allclose(computeNormal(node1Coords, node2Coords, node3Coords), zmaxNormal):
            return True
        return False

    def isboundary(node1Coords, node2Coords, node3Coords, rve):
        """Determines whether a triangle (defined by its 3 nodes) is on the boundary on a rve"""
        if isxmin(node1Coords, node2Coords, node3Coords, rve) or isxmax(node1Coords, node2Coords, node3Coords, rve) or isymin(node1Coords, node2Coords, node3Coords, rve) or isymax(node1Coords, node2Coords, node3Coords, rve) or iszmin(node1Coords, node2Coords, node3Coords, rve) or iszmax(node1Coords, node2Coords, node3Coords, rve):
            return True
        return False

    boundaryElementsTagsToKeep = []
    for i in range(len(boundaryElementsTags[0])):
        if isboundary(boundaryElementsNodesCoords[3*i], boundaryElementsNodesCoords[3*i + 1], boundaryElementsNodesCoords[3*i + 2], rve):
            boundaryElementsTagsToKeep.append(int(boundaryElementsTags[0][i] - ntets))

    nTagsToKeep = len(boundaryElementsTagsToKeep)

    # write files

    gmsh.write(mesh_file[:-5] + "_triangles.mesh")
    gmsh.finalize()

    with open(mesh_file[:-5] + "_triangles.mesh", "r") as file:
        lines = file.readlines()

    with open(output_file, "w+") as outfile:
        outfile.writelines(lines[:-1])
        outfile.write('RequiredTriangles\n')
        outfile.write(str(nTagsToKeep) + '\n')
        for tag in boundaryElementsTagsToKeep:
            outfile.write(str(tag) + '\n')
        outfile.write('End\n')

