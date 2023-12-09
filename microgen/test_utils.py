import numpy as np

def is_periodic(nodes_coords: np.ndarray, tol: float = 1e-8, dim: int = 3) -> bool:

    # bounding box
    xmax = np.max(nodes_coords[:, 0])
    xmin = np.min(nodes_coords[:, 0])
    ymax = np.max(nodes_coords[:, 1])
    ymin = np.min(nodes_coords[:, 1])
    if dim == 3:
        zmax = np.max(nodes_coords[:, 2])
        zmin = np.min(nodes_coords[:, 2])

    # extract face nodes
    face_Xm = np.where(np.abs(nodes_coords[:, 0] - xmin) < tol)[0]
    face_Xp = np.where(np.abs(nodes_coords[:, 0] - xmax) < tol)[0]

    if dim > 1:
        face_Ym = np.where(np.abs(nodes_coords[:, 1] - ymin) < tol)[0]
        face_Yp = np.where(np.abs(nodes_coords[:, 1] - ymax) < tol)[0]

    if dim > 2:  # or dim == 3
        face_Zm = np.where(np.abs(nodes_coords[:, 2] - zmin) < tol)[0]
        face_Zp = np.where(np.abs(nodes_coords[:, 2] - zmax) < tol)[0]

        # sort adjacent faces to ensure node correspondence
    if nodes_coords.shape[1] == 2:  # 2D mesh
        face_Xm = face_Xm[np.argsort(nodes_coords[face_Xm, 1])]
        face_Xp = face_Xp[np.argsort(nodes_coords[face_Xp, 1])]
        if dim > 1:
            face_Ym = face_Ym[np.argsort(nodes_coords[face_Ym, 0])]
            face_Yp = face_Yp[np.argsort(nodes_coords[face_Yp, 0])]

    elif nodes_coords.shape[1] > 2:
        decimal_round = int(-np.log10(tol) - 1)        
        def _sort_dim(indices: np.ndarray, dim_a: int, dim_b: int) -> np.ndarray:
            return indices[
                np.lexsort(
                    (
                        nodes_coords[indices, dim_a],
                        nodes_coords[indices, dim_b].round(decimal_round),
                    )
                )
            ]

        face_Xm = _sort_dim(face_Xm, dim_a=1, dim_b=2)
        face_Xp = _sort_dim(face_Xp, dim_a=1, dim_b=2)
        if dim > 1:
            face_Ym = _sort_dim(face_Ym, dim_a=0, dim_b=2)
            face_Yp = _sort_dim(face_Yp, dim_a=0, dim_b=2)
        if dim > 2:
            face_Zm = _sort_dim(face_Zm, dim_a=0, dim_b=1)
            face_Zp = _sort_dim(face_Zp, dim_a=0, dim_b=1)

    # ==========================
    # test if mesh is periodic:
    # ==========================

    # test if same number of nodes in adjacent faces
    if len(face_Xm) != len(face_Xp):
        return False
    if dim > 1 and len(face_Ym) != len(face_Yp):
        return False
    if dim > 2 and (len(face_Zm) != len(face_Zp)):
        return False

    # check nodes position
    if (nodes_coords[face_Xp, 1:] - nodes_coords[face_Xm, 1:] > tol).any():
        return False
    if (
        dim > 1
        and (nodes_coords[face_Yp, ::2] - nodes_coords[face_Ym, ::2] > tol).any()
    ):
        return False
    if dim > 2 and (nodes_coords[face_Zp, :2] - nodes_coords[face_Zm, :2] > tol).any():
        return False

    return True
