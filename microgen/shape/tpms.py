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
import progressbar

from microgen.operations import fuseShapes
from microgen.shape.basicGeometry import BasicGeometry

SHEET = "sheet"
SKELETAL = "skeletal"
ALL = "all"

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
        surface_function: Callable = None,
        offset: Union[float, Callable] = 0.,
        cell_size: Union[float, np.ndarray] = np.array([1., 1., 1.]),
        repeat_cell: Union[int, np.ndarray] = np.array([1, 1, 1]),
    ) -> None:
        """
        :param center: center of the geometry
        :param orientation: orientation of the geometry
        :param surface_function: tpms function or custom function (f(x, y, z, t) = 0)
        :param type_part: 'sheet' or 'skeletal'
        :param offset: offset of the isosurface to generate thickness
        :param cell_size: By default, the tpms is generated for (1, 1, 1) dimensions but this can be modified by passing 'cell_size' scaling parameter (float or list of float for each dimension)
        :param repeat_cell: By default, the tpms is generated for one unit_cell. 'repeat_cell' parameter allows to repeat the geometry in the three dimensions
        """
        super().__init__(shape="TPMS", center=center, orientation=orientation)

        self.surface_function = surface_function
        self.offset = offset

        if type(cell_size) in [float, int]:
            self.cell_size = np.array([cell_size, cell_size, cell_size])
        else:
            self.cell_size = np.array(cell_size)

        if type(repeat_cell) == int:
            self.repeat_cell = np.array([repeat_cell, repeat_cell, repeat_cell])
        else:
            self.repeat_cell = np.array(repeat_cell)

    def _computeSurfaceFunction(
        self,
        surface_function: Callable,
        resolution: int = 20
    ) -> tuple[pv.UniformGrid, np.ndarray]:

        n_points = resolution * self.repeat_cell 
        grid = pv.UniformGrid()
        grid.dimensions = n_points
        grid.origin = -0.5 * self.cell_size
        grid.spacing = self.cell_size / (n_points - 1)
        
        x, y, z = (grid.points * 2 * np.pi * self.repeat_cell / self.cell_size).T

        surface = surface_function(x, y, z)

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

    def _computeSurfaceVtk(
        self,
        grid : pv.UniformGrid,
        surface: np.ndarray,
        smoothing: int = 100
    ) -> pv.PolyData:
        # mesh = grid.contour(
        #     1, scalars=surface, method="flying_edges", rng=(isovalue, 1)
        # )
        mesh = grid.contour(isosurfaces=[0.], scalars=surface)
        if smoothing > 0:
            mesh = mesh.smooth(n_iter=smoothing)
        mesh.clean(inplace=True)

        return mesh

    def _createShell(
        self,
        mesh: pv.PolyData,
        verbose: bool
    ) -> cq.Shell:
        list_of_triangles = mesh.faces.reshape(-1, 4)[:, 1:]
        list_of_triangles = np.c_[list_of_triangles, list_of_triangles[:, 0]]

        pbar = progressbar.ProgressBar(max_value=len(list_of_triangles))
        if not verbose: pbar = progressbar.NullBar

        faces = []
        for i, ixs in enumerate(list_of_triangles):
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
            pbar.update(i)
        
        if verbose: print("\nGenerating shell surface\n")
        return cq.Shell.makeShell(faces)

    def _createSurface(
        self,
        grid: pv.UniformGrid, 
        surface: np.ndarray, 
        smoothing: int, 
        verbose: bool
    ) -> cq.Shell:
        """
        Create TPMS surface for the corresponding isovalue, return a cq.Shell

        :param isovalue: height isovalue of the given tpms function
        :param resolution: grid resolution for marching cubes algorithm
        :param smoothing: smoothing loop iterations
        """
        mesh = self._computeSurfaceVtk(grid=grid,
                                       surface=surface,
                                       smoothing=smoothing)
        return self._createShell(mesh=mesh, verbose=verbose)

    def _createSurfaces(
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
        grid, surface = self._computeSurfaceFunction(self.surface_function, resolution)

        shells = []
        for isovalue in isovalues:
            mesh = self._computeSurfaceVtk(grid, surface, isovalue, smoothing)
            shells.append(self._createShell(mesh))
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
        grid, surface = self._computeSurfaceFunction(self.surface_function, resolution)
        
        return self._computeSurfaceVtk(grid, surface, isovalue, smoothing)

    def generateSurface(
        self,
        isovalue: float = 0.0,
        resolution: int = 20,
        smoothing: int = 100,
    ) -> cq.Shape:

        shell = self._createSurface(isovalue=isovalue, resolution=resolution, smoothing=smoothing)

        return cq.Shape(shell.wrapped)

    def _get_sheet(self, listShapes):
        sheet = [shape for (n_shapes, shape) in listShapes if n_shapes > 1]
        # to_fuse = [cq.Shape(shape.wrapped) for shape in sheet]
        return fuseShapes(cqShapeList=sheet, retain_edges=True)

    def _get_skeletals(self, listShapes):
        skeletals = [shape for (n_shapes, shape) in listShapes if n_shapes == 1]
        return fuseShapes(cqShapeList=skeletals, retain_edges=False)


    def _separate_sheet_and_skeletal_parts(self,
                                           type_part: str,
                                           solids: list[cq.Workplane],
                                           grid: pv.UniformGrid,
                                           surface_function: np.ndarray,
                                           isosurface: Union[float, np.ndarray],
                                           smoothing: int,
                                           verbose: bool):
        if verbose: print(f"Extracting {type_part} part(s) from generated solids\n")

        if verbose: print("Converting upper test surface")
        upper_test_surface_function = surface_function + isosurface / 6.
        upper_test_surface = self._createSurface(grid=grid,
                                                 surface=upper_test_surface_function, 
                                                 smoothing=smoothing, 
                                                 verbose=verbose)

        if verbose: print("Converting lower test surface")
        lower_test_surface_function = surface_function - isosurface / 6.
        lower_test_surface = self._createSurface(grid=grid,
                                                 surface=lower_test_surface_function, 
                                                 smoothing=smoothing, 
                                                 verbose=verbose)

        
        if verbose: print("Splitting solids\n")
        listShapes: list[tuple[int, cq.Shape]] = []
        for solid in solids:
            # if solid is splitted in several parts by test surfaces then the solid is in the sheet part
            shape = solid.val()
            n_shapes = (solid.split(upper_test_surface)
                             .split(lower_test_surface)
                             .solids().size())
            listShapes.append((n_shapes, shape))

        if verbose: print("Return required part(s)")
        if type_part == SHEET:
            shape = self._get_sheet(listShapes)
        elif type_part == SKELETAL:
            shape = self._get_skeletals(listShapes)
        elif type_part == ALL:
            shape = (self._get_sheet(listShapes), self._get_skeletals(listShapes))
        return shape

    def generate(
        self,
        type_part: str,
        resolution: int = 20,
        smoothing: int = 100,
        verbose: bool=True,
        preview: bool=False
    ) -> cq.Shape:
        """
        Creates thick TPMS geometry (sheet or skeletal part) from surface

        :param resolution: surface file name
        :param smoothing: smoothing loop iterations
        """
        
        if type_part not in [SHEET, SKELETAL, ALL]:
            raise ValueError("type_part must be 'sheet' or 'skeletal' or 'all'")

        n_points = resolution * self.repeat_cell 
        grid = pv.UniformGrid()
        grid.dimensions = n_points
        grid.origin = -0.5 * self.cell_size
        grid.spacing = self.cell_size / (n_points - 1)
        
        x, y, z = (grid.points * 2 * np.pi * self.repeat_cell / self.cell_size).T
        
        surface_function = self.surface_function(x, y, z)
        # f(x, y, z) = isosurface
        isosurface = 2 * np.pi
        if isinstance(self.offset, Callable):
            isosurface *= self.offset(x, y, z)
        else:
            isosurface *= self.offset

        upper_surface_function = surface_function + isosurface / 2.
        mesh_upper_surface = self._computeSurfaceVtk(grid=grid,
                                                     surface=upper_surface_function,
                                                     smoothing=smoothing)
        
        lower_surface_function = surface_function - isosurface / 2.
        mesh_lower_surface = self._computeSurfaceVtk(grid=grid,
                                                     surface=lower_surface_function,
                                                     smoothing=smoothing)

        if preview: 
            pl = pv.Plotter()
            pl.add_mesh(mesh_upper_surface, color="w")
            pl.add_mesh(mesh_lower_surface, color="w")
            pl.add_axes()
            pl.view_xz()
            pl.show()
            return
                                           
        if verbose: print("Converting upper surface")          
        upper_surface = self._createShell(mesh=mesh_upper_surface, verbose=verbose)
        
        if verbose: print("Converting lower surface")
        lower_surface = self._createShell(mesh=mesh_lower_surface, verbose=verbose)

        
        if verbose: print("Splitting box with generated surfaces\n")
        box = cq.Workplane("front").box(*self.cell_size)
        splitted = box.split(upper_surface).split(lower_surface)
        solids: list[cq.Workplane] = splitted.solids().all()

        return self._separate_sheet_and_skeletal_parts(type_part,
                                                       solids,
                                                       grid,
                                                       surface_function,
                                                       isosurface,
                                                       smoothing=smoothing,
                                                       verbose=verbose)


    def generateVtk(
        self,
        resolution: int = 20,
        smoothing: int = 100,
        verbose: bool=True
    ) -> pv.PolyData:
        """
        Creates thick TPMS geometry (sheet or skeletal part) from surface
        Calls generate function and converts cq.Shape to pv.Polydata

        :param resolution: surface file name
        :param smoothing: smoothing loop iterations
        """
        shape = self.generate(resolution=resolution, smoothing=smoothing, verbose=verbose)
        return pv.PolyData(
            shape.toVtkPolyData(tolerance=0.01, angularTolerance=0.1, normals=True)
        )


def gyroid(x: float, y: float, z: float) -> float:
    """
    .. math:: 
       sin(2 \pi x) cos(2 \pi y) + sin(2 \pi y) cos(2 \pi z) + sin(2 \pi z) cos(2 \pi x) = 0

    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.gyroid,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white')
    """
    return (
        sin(x) * cos(y)
        + sin(y) * cos(z)
        + sin(z) * cos(x)
    )


def schwarzP(x: float, y: float, z: float) -> float:
    """
    .. math::
       cos(2 \pi x) + cos(2 \pi y) + cos(2 \pi z) = 0

    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.schwarzP,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white')    
    """
    return cos(x) + cos(y) + cos(z)


def schwarzD(x: float, y: float, z: float) -> float:
    """
    .. math::
       \displaylines{sin(2 \pi x) sin(2 \pi y) sin(2 \pi z) + \\\ sin(2 \pi x) cos(2 \pi y) cos(2 \pi z) + \\\ cos(2 \pi x) sin(2 \pi y) cos(2 \pi z) + \\\ cos(2 \pi x) cos(2 \pi y) sin(2 \pi z) = 0}

    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.schwarzD,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    a = sin(x) * sin(y) * sin(z)
    b = sin(x) * cos(y) * cos(z)
    c = cos(x) * sin(y) * cos(z)
    d = cos(x) * cos(y) * sin(z)
    return a + b + c + d


def neovius(x: float, y: float, z: float) -> float:
    """
    .. math::
       \displaylines{3 cos(2 \pi x) + cos(2 \pi y) + cos(2 \pi z) + \\\ 4 cos(2 \pi x) cos(2 \pi y) cos(2 \pi z) = 0}
           
    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.neovius,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    a = 3 * cos(x) + cos(y) + cos(z)
    b = 4 * cos(x) * cos(y) * cos(z)

    return a + b


def schoenIWP(x: float, y: float, z: float) -> float:
    """
    .. math::
       \displaylines{2 ( cos(2 \pi x) cos(2 \pi y) + \\\ cos(2 \pi y) cos(2 \pi z) + \\\ cos(2 \pi z) cos(2 \pi x)) - \\\ (cos(4 \pi x) + cos(4 \pi y) + cos(4 \pi z)) = 0}
    
    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.schoenIWP,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    a = 2 * (
        cos(x) * cos(y)
        + cos(y) * cos(z)
        + cos(z) * cos(x)
    )
    b = cos(4 * pi * x) + cos(4 * pi * y) + cos(4 * pi * z)

    return a - b


def schoenFRD(x: float, y: float, z: float) -> float:
    """
    .. math::
       \displaylines{4 cos(2 \pi x) cos(2 \pi y) cos(2 \pi z) - \\\ (cos(4 \pi x) cos(4 \pi y) + \\\ cos(4 \pi y) cos(4 \pi z) + \\\ cos(4 \pi z) cos(4 \pi x)) = 0}
            
    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.schoenFRD,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    a = 4 * cos(x) * cos(y) * cos(z)
    b = (
        cos(2 * x) * cos(2 * y)
        + cos(2 * y) * cos(2 * z)
        + cos(2 * z) * cos(2 * x)
    )
    return a - b


def fischerKochS(x: float, y: float, z: float) -> float:
    """
    .. math::
       \displaylines{cos(4 \pi x) sin(2 \pi y) cos(2 \pi z) + \\\ cos(2 \pi x) cos(4 \pi y) sin(2 \pi z) + \\\ sin(2 \pi x) cos(2 \pi y) cos(4 \pi z) = 0}
    
    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.fischerKochS,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    a = cos(2 * x) * sin(y) * cos(z)
    b = cos(x) * cos(2 * y) * sin(z)
    c = sin(x) * cos(y) * cos(2 * z)

    return a + b + c


def pmy(x: float, y: float, z: float) -> float:
    """
    .. math::
       \displaylines{2 cos(2 \pi x) cos(2 \pi y) cos(2 \pi z) + \\\ sin(4 \pi x) sin(2 \pi y) + \\\ sin(2 \pi x) sin(4 \pi z) + \\\ sin(4 \pi y) sin(2 \pi z) = 0}
           
    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.pmy,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    a = 2 * cos(x) * cos(y) * cos(z)
    b = sin(2 * x) * sin(y)
    c = sin(x) * sin(2 * z)
    d = sin(2 * y) * sin(z)

    return a + b + c + d


def honeycomb(x: float, y: float, z: float) -> float:
    """
    .. math::
       sin(2 \pi x) cos(2 \pi y) + sin(2 \pi y) + cos(2 \pi z) = 0
       
    .. jupyter-execute::
       :hide-code:
   
       import microgen

       geometry = microgen.Tpms(
           surface_function=microgen.tpms.honeycomb,
           type_part="sheet"
       )
       shape = geometry.generateSurfaceVtk()

       shape.plot(color='white') 
    """
    return sin(x) * cos(y) + sin(y) + cos(z)


def lidinoid(x: float, y: float, z: float) -> float:
    return 0.5 * (sin(2 * x) * cos(y) * sin(z) + 
                  sin(2 * y) * cos(z) * sin(x) + 
                  sin(2 * z) * cos(x) * sin(y)) - \
           0.5 * (cos(2 * x) * cos(2 * y) +
                  cos(2 * y) * cos(2 * z) + 
                  cos(2 * z) * cos(2 * x)) + 0.3


def split_p(x: float, y: float, z: float) -> float:
    return 1.1 * (sin(2 * x) * cos(y) * sin(z) + 
                  sin(2 * y) * cos(z) * sin(x) + 
                  sin(2 * z) * cos(x) * sin(y)) - \
           0.2 * (cos(2 * x) * cos(2 * y) +
                  cos(2 * y) * cos(2 * z) + 
                  cos(2 * z) * cos(2 * x)) - \
           0.4 * (cos(2 * x) + cos(2 * y) + cos(2 * z))