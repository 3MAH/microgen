from microgen.Functions import fuseParts
from microgen.TriplyPeriodicSurfaces import generateTPMS, is_function_valid
import cadquery as cq

from OCP.StlAPI import StlAPI_Reader

from OCP.TopoDS import TopoDS_Shape


# ----------TPMS-----------------------------------------------------------------------------------------#


class Tpms:
    def __init__(self, center, angle, type_surface, type_part, thickness, number, function=""):
        self.center = center
        self.angle = angle
        self.number = number
        self.name_part = "tpms" + str(self.number)
        self.type_surface = type_surface
        self.type_part = type_part
        self.thickness = thickness
        self.function = function

        if type_surface == "custom":
            if function == "":
                raise ValueError("No given function to evaluate")
            elif not is_function_valid(function):
                raise ValueError(function + " function not valid")



    def createSurfaces(
        self, rve, sizeMesh=0.05, minFacetAngle=10, maxRadius=0.05, path_data=""
    ):
        generateTPMS(
            self.type_surface,
            self.thickness,
            rve,
            sizeMesh,
            minFacetAngle,
            maxRadius,
            path_data,
            function=self.function
        )

    def createTpms(self, path_data, rve):

        if rve is None:
            print("Please add an RVE to generate the TPMS")
        surf_tp = TopoDS_Shape()
        surf_tm = TopoDS_Shape()
        surf_p = TopoDS_Shape()
        surf_m = TopoDS_Shape()
        stl_reader = StlAPI_Reader()
        if path_data != "":
            stl_reader.Read(surf_tp, path_data + "/" + "tpms_testplus.stl")
            stl_reader.Read(surf_tm, path_data + "/" + "tpms_testminus.stl")
            stl_reader.Read(surf_p, path_data + "/" + "tpms_plus.stl")
            stl_reader.Read(surf_m, path_data + "/" + "tpms_minus.stl")
        else:
            stl_reader.Read(surf_tp, "tpms_testplus.stl")
            stl_reader.Read(surf_tm, "tpms_testminus.stl")
            stl_reader.Read(surf_p, "tpms_plus.stl")
            stl_reader.Read(surf_m, "tpms_minus.stl")

        face_cut_tp = cq.Face(surf_tp)
        face_cut_tm = cq.Face(surf_tm)
        face_cut_p = cq.Face(surf_p)
        face_cut_m = cq.Face(surf_m)

        box = cq.Workplane("front").box(rve.dx, rve.dy, rve.dz)

        boxCut = box.split(face_cut_p)
        boxCut = boxCut.split(face_cut_m)

        boxSolids = boxCut.solids().all()
        boxSolidsSize = boxCut.solids().size()

        print("boxSolids", boxSolids)
        print("boxSolidsSize", boxSolidsSize)
        print("boxCut", boxCut)

        listSolids = []

        for solid in boxSolids:
            temp = solid.split(face_cut_tp)
            temp = temp.split(face_cut_tm)
            tempSize = temp.solids().size()
            listSolids.append((tempSize, solid.val()))

        sheet = [el[1] for el in listSolids if el[0] > 1]
        skeletal = [el[1] for el in listSolids if el[0] == 1]

        print("sheet", sheet)
        print("skeletal", skeletal)

        if self.type_part == "sheet":
            to_fuse = [cq.Shape(s.wrapped) for s in sheet]
            return_object = fuseParts(to_fuse, True)
            return cq.Workplane().add(return_object[0])
        elif self.type_part == "skeletal":
            to_fuse = [cq.Shape(s.wrapped) for s in skeletal]
            return_object = fuseParts(to_fuse, False)
            return cq.Workplane().add(return_object[0])
