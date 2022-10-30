"""
=============================================
TPMS (:mod:`microgen.shape.tpms`)
=============================================

.. jupyter-execute::
   :hide-code:

   pyvista.global_theme.smooth_shading = True

"""
from typing import Callable, Union

import numpy as np
import cadquery as cq
import pyvista as pv
from numpy import cos, pi, sin

from microgen.operations import fuseShapes, rescale, repeatShape, repeatPolyData
from microgen.rve import Rve
from microgen.shape.basicGeometry import BasicGeometry

SHEET = "sheet"
SKELETAL = "skeletal"

class Tpms(BasicGeometry):
    """
    Class to generate Triply Periodical Minimal Surfaces (TPMS)
    geometry from a given mathematical function, with given thickness

    functions available :
        - :class:`~microgen.shape.tpms.gyroid`
        - :class:`~microgen.shape.tpms.schwarzP`
        - :class:`~microgen.shape.tpms.schwarzD`
        - :class:`~microgen.shape.tpms.neovius`
        - :class:`~microgen.shape.tpms.schoenIWP`
        - :class:`~microgen.shape.tpms.schoenFRD`
        - :class:`~microgen.shape.tpms.fischerKochS`
        - :class:`~microgen.shape.tpms.pmy`
        - :class:`~microgen.shape.tpms.honeycomb`
    """

    def __init__(
        self,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        surface_function: Callable[[float, float, float, float], float] = None,
        type_part: str = SHEET,
        thickness: float = 0,
        cell_size: Union[float, tuple[float, float, float]] = (1., 1., 1.),
        repeat_cell: Union[int, tuple[int, int, int]] = (1, 1, 1),
    ) -> None:
        """
        :param center: center of the geometry
        :param orientation: orientation of the geometry
        :param surface_function: tpms function or custom function (f(x, y, z, t) = 0)
        :param type_part: 'sheet' or 'skeletal'
        :param thickness: thickness of the tpms
        :param cell_size: By default, the tpms is generated for (1, 1, 1) dimensions but this can be modified by passing 'cell_size' scaling parameter (float or list of float for each dimension)
        :param repeat_cell: By default, the tpms is generated for one unit_cell. 'repeat_cell' parameter allows to repeat the geometry in the three dimensions
        """
        super().__init__(shape="TPMS", center=center, orientation=orientation)

        self.surface_function = surface_function

        if type_part not in [SHEET, SKELETAL]:
            raise ValueError("type_part must be 'sheet' or 'skeletal'")
        self.type_part = type_part

        self.thickness = thickness * 2 * np.pi

        if type(cell_size) in [float, int]:
            self.cell_size = (cell_size, cell_size, cell_size)
        else:
            self.cell_size = cell_size

        if type(repeat_cell) == int:
            self.repeat_cell = (repeat_cell, repeat_cell, repeat_cell)
        else:
            self.repeat_cell = repeat_cell

    def computeSurfaceFunction(
        self,
        surface_function: Callable,
        resolution: int = 20
    ) -> tuple[pv.UniformGrid, np.ndarray]:

        x_min, y_min, z_min = -0.5 * self.cell_size[0], -0.5 * self.cell_size[1], -0.5 * self.cell_size[2]
        grid = pv.UniformGrid(
            dims=(resolution, resolution, resolution),
            spacing=(self.cell_size[0] / (resolution - 1),
                     self.cell_size[1] / (resolution - 1),
                     self.cell_size[2] / (resolution - 1)),
            origin=(x_min, y_min, z_min),
        )
        x, y, z = grid.points.T

        surface = surface_function(x, y, z, self.repeat_cell, self.cell_size)

        return grid, surface

    # def computeSurfaceFunctionOnCylinder(
    #     self,
    #     resolution: int
    # ) -> tuple[pv.StructuredGrid, np.ndarray]:
    #     res = 100j
    #     a, b, c = lattice_params = 1, 1, 1
    #     kx, ky, kz = [2*np.pi/lattice_param for lattice_param in lattice_params]
    #     r_aux, phi, z = np.mgrid[0:a:res, 0:b:res, 0:3*c:res]

    #     # convert r_aux range to actual radii
    #     r1, r2 = 1.5, 2
    #     r = r2/a*r_aux + r1/a*(1 - r_aux)

    #     surface = self.surface_function(r_aux, phi / 12, z, self.repeat_cell, self.cell_size)

    #     x = r * np.cos(phi*kx)
    #     y = r * np.sin(phi*ky)
    #     grid = pv.StructuredGrid(x, y, z)

    #     return grid, surface

    # def cylinder(self, isovalue, resolution, smoothing):
    #     grid, surface = self.computeSurfaceFunctionOnCylinder(resolution)
        
    #     return self.computeSurfaceVtk(grid, surface, isovalue, smoothing)

    def computeSurfaceVtk(
        self,
        grid : pv.UniformGrid,
        surface,
        isovalue: float = 0,
        smoothing: int = 100
    ) -> pv.PolyData:
        # mesh = grid.contour(
        #     1, scalars=surface, method="flying_edges", rng=(isovalue, 1)
        # )
        mesh = grid.contour(isosurfaces=[isovalue], scalars=surface)
        if smoothing > 0:
            mesh = mesh.smooth(n_iter=smoothing)
        mesh.clean(inplace=True)

        return mesh

    def createShell(
        self,
        mesh: pv.PolyData
    ) -> cq.Shell:
        list_of_triangles = mesh.faces.reshape(-1, 4)[:, 1:]
        list_of_triangles = np.c_[list_of_triangles, list_of_triangles[:, 0]]

        faces = []
        for ixs in list_of_triangles:
            lines = []
            for v1, v2 in zip(ixs[:], ixs[1:]):
                vertice_coords1 = mesh.points[v1]
                vertice_coords2 = mesh.points[v2]
                lines.append(
                    cq.Edge.makeLine(
                        cq.Vector(*vertice_coords1), cq.Vector(*vertice_coords2)
                    )
                )
            wire = cq.Wire.assembleEdges(lines)
            faces.append(cq.Face.makeFromWires(wire))
        return cq.Shell.makeShell(faces)

    def createSurface(
        self,
        isovalue: float = 0,
        resolution: int = 20,
        smoothing: int = 100,
    ) -> cq.Shell:
        """
        Create TPMS surface for the corresponding isovalue, return a cq.Shell

        :param isovalue: height isovalue of the given tpms function
        :param resolution: grid resolution for marching cubes algorithm
        :param smoothing: smoothing loop iterations
        """
        grid, surface = self.computeSurfaceFunction(self.surface_function, resolution)

        mesh = self.computeSurfaceVtk(grid, surface, isovalue, smoothing)

        return self.createShell(mesh)

    def createSurfaces(
        self,
        isovalues: list[float],
        resolution: int = 20,
        smoothing: int = 100,
    ) -> list[cq.Shell]:
        """
        Create TPMS surfaces for the corresponding isovalue, return a list of cq.Shell

        :param numsber_surfaces: number of surfaces
        :param isovalues: height isovalues of the given tpms function
        :param resolution: surface file name
        :param smoothing: smoothing loop iterations
        """
        grid, surface = self.computeSurfaceFunction(self.surface_function, resolution)

        shells = []
        for isovalue in isovalues:
            mesh = self.computeSurfaceVtk(grid, surface, isovalue, smoothing)
            shells.append(self.createShell(mesh))

        return shells

    def generateSurfaceVtk(
        self,
        isovalue: float = 0,
        resolution: int = 20,
        smoothing: int = 100,
    ) -> pv.PolyData:
        """
        Create TPMS surface for the corresponding isovalue, returns a pv.Polydata

        example:
            shape = geometry.generateSurfaceVtk(isovalue=0, resolution=100, smoothing=100)
            shape.plot()

        :param isovalue: height isovalue of the given tpms function
        :param resolution: surface file name
        :param smoothing: smoothing loop iterations
        """
        grid, surface = self.computeSurfaceFunction(self.surface_function, resolution)
        
        return self.computeSurfaceVtk(grid, surface, isovalue, smoothing)

    def generateSurface(
        self,
        isovalue: float = 0.0,
        resolution: int = 20,
        smoothing: int = 100,
    ) -> cq.Shape:

        shell = self.createSurface(isovalue=isovalue, resolution=resolution, smoothing=smoothing)

        return cq.Shape(shell.wrapped)

    def generate(
        self,
        resolution: int = 20,
        smoothing: int = 100,
    ) -> cq.Shape:
        """
        Creates thick TPMS geometry (sheet or skeletal part) from surface

        :param resolution: surface file name
        :param smoothing: smoothing loop iterations
        """

        isovalues = [
            -self.thickness / 2.,
            -self.thickness / 6.0,
            self.thickness / 6.0,
            self.thickness / 2.,
        ]
        shells = self.createSurfaces(isovalues=isovalues, resolution=resolution, smoothing=smoothing)

        face_cut_tp = shells[2]
        face_cut_tm = shells[1]
        face_cut_p = shells[3]
        face_cut_m = shells[0]

        box_wp = cq.Workplane("front").box(1, 1, 1)

        boxCut_wp = box_wp.split(face_cut_p)
        boxCut_wp = boxCut_wp.split(face_cut_m)

        boxWorkplanes: list[cq.Workplane] = boxCut_wp.solids().all()

        listShapes: list[tuple[int, cq.Shape]] = []
        for wp in boxWorkplanes:
            n_shapes = wp.split(face_cut_tp).split(face_cut_tm).solids().size()
            shape = wp.val()
            listShapes.append((n_shapes, shape))

        if self.type_part == "sheet":
            sheet = [shape for (n_shapes, shape) in listShapes if n_shapes > 1]
            to_fuse = [cq.Shape(shape.wrapped) for shape in sheet]
            shape = fuseShapes(to_fuse, True)
        elif self.type_part == "skeletal":
            skeletal = [shape for (n_shapes, shape) in listShapes if n_shapes == 1]
            shape = skeletal[0]
            # to_fuse = [cq.Shape(shape.wrapped) for shape in skeletal]
            # shape = fuseShapes(to_fuse, False)
        return shape

    def generateVtk(
        self,
        resolution: int = 20,
        smoothing: int = 100,
    ) -> pv.PolyData:
        """
        Creates thick TPMS geometry (sheet or skeletal part) from surface
        Calls generate function and converts cq.Shape to pv.Polydata

        :param resolution: surface file name
        :param smoothing: smoothing loop iterations
        """
        shape = self.generate(resolution=resolution, smoothing=smoothing)
        return pv.PolyData(
            shape.toVtkPolyData(tolerance=0.01, angularTolerance=0.1, normals=True)
        )


def lidinoid(x: float, y: float, z: float, q: tuple[int, int, int], l: tuple[float, float, float]) -> float:
    kx = 2* np.pi * q[0] / l[0]
    ky = 2* np.pi * q[1] / l[1]
    kz = 2* np.pi * q[2] / l[2]
    return 0.5*(sin(2*kx*x)*cos(ky*y)*sin(kz*z) + 
                sin(2*ky*y)*cos(kz*z)*sin(kx*x) + 
                sin(2*kz*z)*cos(kx*x)*sin(ky*y)) - \
           0.5*(cos(2*kx*x)*cos(2*ky*y) +
                cos(2*ky*y)*cos(2*kz*z) + 
                cos(2*kz*z)*cos(2*kx*x)) + 0.15


def gyroid(x: float, y: float, z: float, q: tuple[int, int, int], l: tuple[float, float, float]) -> float:
    """
    .. math:: 
       sin(2 \pi x) cos(2 \pi y) + sin(2 \pi y) cos(2 \pi z) + sin(2 \pi z) cos(2 \pi x) = 0

    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.gyroid,
           type_part="sheet",
           thickness=0.05
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white')
    """
    kx = 2* np.pi * q[0] / l[0]
    ky = 2* np.pi * q[1] / l[1]
    kz = 2* np.pi * q[2] / l[2]
    return (
        sin(kx * x) * cos(ky * y)
        + sin(ky * y) * cos(kz * z)
        + sin(kz * z) * cos(kx * x)
    )


def schwarzP(x: float, y: float, z: float, q: tuple[int, int, int], l: tuple[float, float, float]) -> float:
    """
    .. math::
       cos(2 \pi x) + cos(2 \pi y) + cos(2 \pi z) = 0

    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.schwarzP,
           type_part="sheet",
           thickness=0.05
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white')    
    """
    kx = 2* np.pi * q[0] / l[0]
    ky = 2* np.pi * q[1] / l[1]
    kz = 2* np.pi * q[2] / l[2]
    return cos(kx * x) + cos(ky * y) + cos(kz * z)


def schwarzD(x: float, y: float, z: float, q: tuple[int, int, int], l: tuple[float, float, float]) -> float:
    """
    .. math::
       \displaylines{sin(2 \pi x) sin(2 \pi y) sin(2 \pi z) + \\\ sin(2 \pi x) cos(2 \pi y) cos(2 \pi z) + \\\ cos(2 \pi x) sin(2 \pi y) cos(2 \pi z) + \\\ cos(2 \pi x) cos(2 \pi y) sin(2 \pi z) = 0}

    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.schwarzD,
           type_part="sheet",
           thickness=0.05
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    kx = 2* np.pi * q[0] / l[0]
    ky = 2* np.pi * q[1] / l[1]
    kz = 2* np.pi * q[2] / l[2]
    a = sin(kx * x) * sin(ky * y) * sin(kz * z)
    b = sin(kx * x) * cos(ky * y) * cos(kz * z)
    c = cos(kx * x) * sin(ky * y) * cos(kz * z)
    d = cos(kx * x) * cos(ky * y) * sin(kz * z)
    return a + b + c + d


def neovius(x: float, y: float, z: float, q: tuple[int, int, int], l: tuple[float, float, float]) -> float:
    """
    .. math::
       \displaylines{3 cos(2 \pi x) + cos(2 \pi y) + cos(2 \pi z) + \\\ 4 cos(2 \pi x) cos(2 \pi y) cos(2 \pi z) = 0}
           
    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.neovius,
           type_part="sheet",
           thickness=0.05
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    kx = 2* np.pi * q[0] / l[0]
    ky = 2* np.pi * q[1] / l[1]
    kz = 2* np.pi * q[2] / l[2]
    a = 3 * cos(kx * x) + cos(ky * y) + cos(kz * z)
    b = 4 * cos(kx * x) * cos(ky * y) * cos(kz * z)

    return a + b


def schoenIWP(x: float, y: float, z: float, q: tuple[int, int, int], l: tuple[float, float, float]) -> float:
    """
    .. math::
       \displaylines{2 ( cos(2 \pi x) cos(2 \pi y) + \\\ cos(2 \pi y) cos(2 \pi z) + \\\ cos(2 \pi z) cos(2 \pi x)) - \\\ (cos(4 \pi x) + cos(4 \pi y) + cos(4 \pi z)) = 0}
    
    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.schoenIWP,
           type_part="sheet",
           thickness=0.05
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    kx = 2* np.pi * q[0] / l[0]
    ky = 2* np.pi * q[1] / l[1]
    kz = 2* np.pi * q[2] / l[2]
    a = 2 * (
        cos(kx * x) * cos(ky * y)
        + cos(ky * y) * cos(kz * z)
        + cos(kz * z) * cos(kx * x)
    )
    b = cos(4 * pi * x) + cos(4 * pi * y) + cos(4 * pi * z)

    return a - b


def schoenFRD(x: float, y: float, z: float, q: tuple[int, int, int], l: tuple[float, float, float]) -> float:
    """
    .. math::
       \displaylines{4 cos(2 \pi x) cos(2 \pi y) cos(2 \pi z) - \\\ (cos(4 \pi x) cos(4 \pi y) + \\\ cos(4 \pi y) cos(4 \pi z) + \\\ cos(4 \pi z) cos(4 \pi x)) = 0}
            
    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.schoenFRD,
           type_part="sheet",
           thickness=0.05
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    kx = 2* np.pi * q[0] / l[0]
    ky = 2* np.pi * q[1] / l[1]
    kz = 2* np.pi * q[2] / l[2]
    a = 4 * cos(kx * x) * cos(ky * y) * cos(kz * z)
    b = (
        cos(2 * kx * x) * cos(2 * ky * y)
        + cos(2 * ky * y) * cos(2 * kz * z)
        + cos(2 * kz * z) * cos(2 * kx * x)
    )
    return a - b


def fischerKochS(x: float, y: float, z: float, q: tuple[int, int, int], l: tuple[float, float, float]) -> float:
    """
    .. math::
       \displaylines{cos(4 \pi x) sin(2 \pi y) cos(2 \pi z) + \\\ cos(2 \pi x) cos(4 \pi y) sin(2 \pi z) + \\\ sin(2 \pi x) cos(2 \pi y) cos(4 \pi z) = 0}
    
    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.fischerKochS,
           type_part="sheet",
           thickness=0.05
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    kx = 2* np.pi * q[0] / l[0]
    ky = 2* np.pi * q[1] / l[1]
    kz = 2* np.pi * q[2] / l[2]
    a = cos(2 * kx * x) * sin(ky * y) * cos(kz * z)
    b = cos(kx * x) * cos(2 * ky * y) * sin(kz * z)
    c = sin(kx * x) * cos(ky * y) * cos(2 * kz * z)

    return a + b + c


def pmy(x: float, y: float, z: float, q: tuple[int, int, int], l: tuple[float, float, float]) -> float:
    """
    .. math::
       \displaylines{2 cos(2 \pi x) cos(2 \pi y) cos(2 \pi z) + \\\ sin(4 \pi x) sin(2 \pi y) + \\\ sin(2 \pi x) sin(4 \pi z) + \\\ sin(4 \pi y) sin(2 \pi z) = 0}
           
    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.pmy,
           type_part="sheet",
           thickness=0.05
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    kx = 2* np.pi * q[0] / l[0]
    ky = 2* np.pi * q[1] / l[1]
    kz = 2* np.pi * q[2] / l[2]
    a = 2 * cos(kx * x) * cos(ky * y) * cos(kz * z)
    b = sin(2 * kx * x) * sin(ky * y)
    c = sin(kx * x) * sin(2 * kz * z)
    d = sin(2 * ky * y) * sin(kz * z)

    return a + b + c + d


def honeycomb(x: float, y: float, z: float, q: tuple[int, int, int], l: tuple[float, float, float]) -> float:
    """
    .. math::
       sin(2 \pi x) cos(2 \pi y) + sin(2 \pi y) + cos(2 \pi z) = 0
       
    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.honeycomb,
           type_part="sheet",
           thickness=0.05
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    kx = 2* np.pi * q[0] / l[0]
    ky = 2* np.pi * q[1] / l[1]
    kz = 2* np.pi * q[2] / l[2]
    return sin(kx * x) * cos(ky * y) + sin(ky * y) + cos(kz * z)