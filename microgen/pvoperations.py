import os

import pyvista as pv
import numpy as np

def rotatePvEuler(
    object: pv.PolyData,
    center: np.ndarray,
    psi: float,
    theta: float,
    phi: float,
) -> pv.PolyData:
    """

    Rotates object according to XZX Euler angle convention

    Parameters
    ----------
    object :
        Object to rotate
    center :
        numpy array (x, y, z)
    psi, theta, phi :
        Euler angles

    Returns
    -------
    object_r :
        Rotated object
    """

    u = (np.cos(psi), np.sin(psi), 0.0)
    z2 = (
            np.sin(psi) * np.sin(theta),
            -np.sin(theta) * np.cos(psi),
            np.cos(theta)
    )

    object_r = object.rotate_vector(
        vector=(0,0,1),
        angle=psi,
        point=tuple(center)
    )
    object_r = object.rotate_vector(
        vector=u,
        angle=theta,
        point=tuple(center)
    )
    object_r = object.rotate_vector(
        vector=z2,
        angle=phi,
        point=tuple(center)
    )
    return object_r

def fusePvParts(
    pvShapeList: list[pv.PolyData]
) -> pv.PolyData:
    """

    Parameters
    ----------
    pvShapeList : list of shapes to fuse

    Returns
    -------
    pv.PolyData : fused object
    """
    for pv_PolyData in pvShapeList:
        if pv_PolyData.is_all_triangles == False:
            pv_PolyData.triangulate()
        
    pv_PolyDatas = pvShapeList[0]
    for i in range(1, len(pvShapeList)):
        fuse = pv_PolyDatas.boolean_union(pvShapeList[i])

    return fuse
