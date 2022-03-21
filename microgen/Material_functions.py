# -*- coding: utf8 -*-
from microgen.Functions import removeEmptyLines
from microgen.Material_sections import MatSection


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
