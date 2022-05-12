# -*- coding: utf8 -*-
from microgen.Operations import removeEmptyLines

import numpy as np


class MatSection:
    def __init__(
        self,
        number,
        name,
        umat_name,
        psi_mat,
        theta_mat,
        phi_mat,
        nprops,
        nstatev,
        props,
    ):
        self.number = number  # Number of the section
        self.name = name  # Name of the section
        self.umat_name = umat_name  # Name of the constitutive law
        self.psi_mat = psi_mat  # Name of the Euler angle for rotation around the z axis
        self.theta_mat = (
            theta_mat  # Name of the Euler angle for rotation around the x' axis
        )
        self.phi_mat = (
            phi_mat  # Name of the Euler angle for rotation around the z'' axis
        )
        self.nprops = nprops  # Number of material parameters
        self.nstatev = (
            nstatev  # Number of internal variables for user constitutive laws
        )
        self.props = props  # Table of the material properties

        c1 = np.cos(self.psi_mat * np.pi / 180.0)
        s1 = np.sin(self.psi_mat * np.pi / 180.0)
        c2 = np.cos(self.theta_mat * np.pi / 180.0)
        s2 = np.sin(self.theta_mat * np.pi / 180.0)
        c3 = np.cos(self.phi_mat * np.pi / 180.0)
        s3 = np.sin(self.phi_mat * np.pi / 180.0)
        self.R = np.array(
            [
                [c3 * c2 * c1 - s3 * s1, -c1 * s3 - c3 * c2 * s1, c3 * s2],
                [c3 * s1 + c2 * c1 * s3, c3 * c1 - c2 * s3 * s1, s3 * s2],
                [-c1 * s2, s2 * s1, c2],
            ]
        )

        x = np.array([[1.0], [0.0], [0.0]])
        y = np.array([[0.0], [1.0], [0.0]])
        x2mat = np.ndarray.flatten(self.R.dot(x))
        y2mat = np.ndarray.flatten(self.R.dot(y))

        self.p0 = (0.0, 0.0, 0.0)
        self.p1 = x2mat.tolist()
        self.p2 = y2mat.tolist()


def readSections(path_data, section_file):

    nsections = 0
    sections = []
    path_inputfile = path_data + "/" + section_file

    removeEmptyLines(path_inputfile)
    table_props = []

    try:
        fp = open(path_inputfile)
        for line in enumerate(fp):
            if line != "\n":
                if nsections > 0:
                    row_split = line[1].split()
                    table_props.append(row_split)
                nsections += 1

    finally:
        fp.seek(0)
        fp.close()

    for row_split in table_props:
        nprops = int(row_split[6])
        props = []
        for prop in range(0, nprops):
            props.append(float(row_split[8 + prop]))

        # sec = msection(number_section,name_section,umat_name,psi_mat,theta_mat,phi_mat,len(props),nstatev,props)
        sec = MatSection(
            int(row_split[0]),
            row_split[1],
            row_split[2],
            float(row_split[3]),
            float(row_split[4]),
            float(row_split[5]),
            nprops,
            int(row_split[7]),
            props,
        )
        sections.append(sec)

    return sections
