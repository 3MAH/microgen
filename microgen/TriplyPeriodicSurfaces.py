import os
import numpy as np
import pygalmesh


class Hyperboloid(pygalmesh.DomainBase):
    def __init__(self, max_edge_size_at_feature_edges):
        super().__init__()
        self.z0 = -1.0
        self.z1 = 1.0
        self.waist_radius = 0.5
        self.max_edge_size_at_feature_edges = max_edge_size_at_feature_edges

    def eval(self, x):
        if self.z0 < x[2] and x[2] < self.z1:
            return x[0] ** 2 + x[1] ** 2 - (x[2] ** 2 + self.waist_radius) ** 2
        return 1.0

    def getBoundingSphereSquaredRadius(self):
        z_max = max(abs(self.z0), abs(self.z1))
        r_max = z_max**2 + self.waist_radius
        return r_max * r_max + z_max * z_max

    def getFeatures(self):
        radius0 = self.z0**2 + self.waist_radius
        n0 = int(2 * np.pi * radius0 / self.max_edge_size_at_feature_edges)
        circ0 = [
            [
                radius0 * np.cos((2 * np.pi * k) / n0),
                radius0 * np.sin((2 * np.pi * k) / n0),
                self.z0,
            ]
            for k in range(n0)
        ]
        circ0.append(circ0[0])

        radius1 = self.z1**2 + self.waist_radius
        n1 = int(2 * np.pi * radius1 / self.max_edge_size_at_feature_edges)
        circ1 = [
            [
                radius1 * np.cos((2 * np.pi * k) / n1),
                radius1 * np.sin((2 * np.pi * k) / n1),
                self.z1,
            ]
            for k in range(n1)
        ]
        circ1.append(circ1[0])
        return [circ0, circ1]


class SchwarzP(pygalmesh.DomainBase):
    def __init__(self, h):
        super().__init__()
        self.h = h
        self.z0 = 0.0
        self.z1 = 1.0
        self.waist_radius = 0.5

    def getBoundingSphereSquaredRadius(self):
        # z_max = max(abs(self.z0), abs(self.z1))
        # r_max = z_max**2 + self.waist_radius
        return 1.1

    def eval(self, x):
        x2 = np.cos(x[0] * 2 * np.pi)
        y2 = np.cos(x[1] * 2 * np.pi)
        z2 = np.cos(x[2] * 2 * np.pi)
        return x2 + y2 + z2 + self.h


class SchwarzD(pygalmesh.DomainBase):
    def __init__(self, h):
        super().__init__()
        self.h = h
        self.z0 = 0.0
        self.z1 = 1.0
        self.waist_radius = 0.5

    def getBoundingSphereSquaredRadius(self):
        # z_max = max(abs(self.z0), abs(self.z1))
        # r_max = z_max**2 + self.waist_radius
        return 1.1

    def eval(self, x):
        a = (
            np.sin(x[0] * 2 * np.pi)
            * np.sin(x[1] * 2 * np.pi)
            * np.sin(x[2] * 2 * np.pi)
        )
        b = (
            np.sin(x[0] * 2 * np.pi)
            * np.cos(x[1] * 2 * np.pi)
            * np.cos(x[2] * 2 * np.pi)
        )
        c = (
            np.cos(x[0] * 2 * np.pi)
            * np.sin(x[1] * 2 * np.pi)
            * np.cos(x[2] * 2 * np.pi)
        )
        d = (
            np.cos(x[0] * 2 * np.pi)
            * np.cos(x[1] * 2 * np.pi)
            * np.sin(x[2] * 2 * np.pi)
        )
        return a + b + c + d + self.h


class Neovius(pygalmesh.DomainBase):
    def __init__(self, h):
        super().__init__()
        self.h = h
        self.z0 = 0.0
        self.z1 = 1.0
        self.waist_radius = 0.5

    def getBoundingSphereSquaredRadius(self):
        # z_max = max(abs(self.z0), abs(self.z1))
        # r_max = z_max**2 + self.waist_radius
        return 1.1

    def eval(self, x):
        a = 3.0 * (
            np.cos(x[0] * 2 * np.pi)
            + np.cos(x[1] * 2 * np.pi)
            + np.cos(x[2] * 2 * np.pi)
        )
        b = 4.0 * (
            np.cos(x[0] * 2 * np.pi)
            * np.cos(x[1] * 2 * np.pi)
            * np.cos(x[2] * 2 * np.pi)
        )
        return a + b + self.h


class Gyroid(pygalmesh.DomainBase):
    def __init__(self, h):
        super().__init__()
        self.h = h
        self.z0 = 0.0
        self.z1 = 1.0
        self.waist_radius = 0.5

    def getBoundingSphereSquaredRadius(self):
        # z_max = max(abs(self.z0), abs(self.z1))
        # r_max = z_max**2 + self.waist_radius
        return 1.1

    def eval(self, x):
        x2 = np.sin(x[0] * 2 * np.pi) * np.cos(x[1] * 2 * np.pi)
        y2 = np.sin(x[1] * 2 * np.pi) * np.cos(x[2] * 2 * np.pi)
        z2 = np.sin(x[2] * 2 * np.pi) * np.cos(x[0] * 2 * np.pi)
        if np.abs(x[0]) + np.abs(x[1]) + np.abs(x[2]) > 1.0e-8:
            return x2 + y2 + z2 + self.h
        else:
            return 1.0


def generateTPMS(
    type_tpms,
    thickness,
    rve,
    sizeMesh=0.05,
    minFacetAngle=10.0,
    maxRadius=0.05,
    path_data="",
):

    # cube = pygalmesh.Cuboid([0, 0, 0], [rve.dx, rve.dy, rve.dz])

    if type_tpms == "gyroid":
        s_testplus = Gyroid(thickness / 2.0)
        s_testminus = Gyroid(-thickness / 2.0)
        s_plus = Gyroid(thickness)
        s_minus = Gyroid(-1.0 * thickness)
    elif type_tpms == "schwarzP":
        s_testplus = SchwarzP(thickness / 2.0)
        s_testminus = SchwarzP(-thickness / 2.0)
        s_plus = SchwarzP(thickness)
        s_minus = SchwarzP(-1.0 * thickness)
    elif type_tpms == "schwarzD":
        s_testplus = SchwarzP(thickness / 2.0)
        s_testminus = SchwarzP(-thickness / 2.0)
        s_plus = SchwarzD(thickness)
        s_minus = SchwarzD(-1.0 * thickness)
    elif type_tpms == "neovius":
        s_testplus = Neovius(thickness / 2.0)
        s_testminus = Neovius(-thickness / 2.0)
        s_plus = Neovius(thickness)
        s_minus = Neovius(-1.0 * thickness)
    else:
        print("Error, the tpms is not recognized")
        return False

    mesh_surf_testplus = pygalmesh.generate_surface_mesh(
        s_testplus,
        min_facet_angle=minFacetAngle,
        max_radius_surface_delaunay_ball=maxRadius,
        max_facet_distance=sizeMesh,
    )

    mesh_surf_testminus = pygalmesh.generate_surface_mesh(
        s_testminus,
        min_facet_angle=minFacetAngle,
        max_radius_surface_delaunay_ball=maxRadius,
        max_facet_distance=sizeMesh,
    )

    mesh_surf_plus = pygalmesh.generate_surface_mesh(
        s_plus,
        min_facet_angle=minFacetAngle,
        max_radius_surface_delaunay_ball=maxRadius,
        max_facet_distance=sizeMesh,
    )

    mesh_surf_minus = pygalmesh.generate_surface_mesh(
        s_minus,
        min_facet_angle=minFacetAngle,
        max_radius_surface_delaunay_ball=maxRadius,
        max_facet_distance=sizeMesh,
    )

    if path_data != '':
        if not(os.path.isdir(path_data)): os.mkdir(path_data)
        mesh_surf_testplus.write(path_data + '/' + 'tpms_testplus.stl')
        mesh_surf_testminus.write(path_data + '/' + 'tpms_testminus.stl')
        mesh_surf_plus.write(path_data + '/' + 'tpms_plus.stl')
        mesh_surf_minus.write(path_data + '/' + 'tpms_minus.stl')
    else:
        mesh_surf_testplus.write("tpms_testplus.stl")
        mesh_surf_testminus.write("tpms_testminus.stl")
        mesh_surf_plus.write("tpms_plus.stl")
        mesh_surf_minus.write("tpms_minus.stl")

    return True
