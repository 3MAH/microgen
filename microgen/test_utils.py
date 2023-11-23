import numpy as np

def is_periodic(nodes_coords, tol=1e-8, dim=3) -> bool:
    # bounding box
    xmax = np.max(nodes_coords[:, 0])
    xmin = np.min(nodes_coords[:, 0])
    ymax = np.max(nodes_coords[:, 1])
    ymin = np.min(nodes_coords[:, 1])
    if dim == 3:
        zmax = np.max(nodes_coords[:, 2])
        zmin = np.min(nodes_coords[:, 2])

    # extract face nodes
    left = np.where(np.abs(nodes_coords[:, 0] - xmin) < tol)[0]
    right = np.where(np.abs(nodes_coords[:, 0] - xmax) < tol)[0]

    if dim > 1:
        bottom = np.where(np.abs(nodes_coords[:, 1] - ymin) < tol)[0]
        top = np.where(np.abs(nodes_coords[:, 1] - ymax) < tol)[0]

    if dim > 2:  # or dim == 3
        back = np.where(np.abs(nodes_coords[:, 2] - zmin) < tol)[0]
        front = np.where(np.abs(nodes_coords[:, 2] - zmax) < tol)[0]

        # sort adjacent faces to ensure node correspondance
    if nodes_coords.shape[1] == 2:  # 2D mesh
        left = left[np.argsort(nodes_coords[left, 1])]
        right = right[np.argsort(nodes_coords[right, 1])]
        if dim > 1:
            bottom = bottom[np.argsort(nodes_coords[bottom, 0])]
            top = top[np.argsort(nodes_coords[top, 0])]

    elif nodes_coords.shape[1] > 2:
        decimal_round = int(-np.log10(tol) - 1)
        left = left[
            np.lexsort(
                (nodes_coords[left, 1], nodes_coords[left, 2].round(decimal_round))
            )
        ]
        right = right[
            np.lexsort(
                (nodes_coords[right, 1], nodes_coords[right, 2].round(decimal_round))
            )
        ]
        if dim > 1:
            bottom = bottom[
                np.lexsort(
                    (
                        nodes_coords[bottom, 0],
                        nodes_coords[bottom, 2].round(decimal_round),
                    )
                )
            ]
            top = top[
                np.lexsort(
                    (nodes_coords[top, 0], nodes_coords[top, 2].round(decimal_round))
                )
            ]
        if dim > 2:
            back = back[
                np.lexsort(
                    (nodes_coords[back, 0], nodes_coords[back, 1].round(decimal_round))
                )
            ]
            front = front[
                np.lexsort(
                    (
                        nodes_coords[front, 0],
                        nodes_coords[front, 1].round(decimal_round),
                    )
                )
            ]

    # ==========================
    # test if mesh is periodic:
    # ==========================

    # test if same number of nodes in adjacent faces
    if len(left) != len(right):
        return False
    if dim > 1 and len(bottom) != len(top):
        return False
    if dim > 2 and (len(back) != len(front)):
        return False

    # check nodes position
    if (nodes_coords[right, 1:] - nodes_coords[left, 1:] > tol).any():
        return False
    if dim > 1 and (nodes_coords[top, ::2] - nodes_coords[bottom, ::2] > tol).any():
        return False
    if dim > 2 and (nodes_coords[front, :2] - nodes_coords[back, :2] > tol).any():
        return False

    return True
