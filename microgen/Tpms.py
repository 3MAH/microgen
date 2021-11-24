from microgen.Functions import *
from microgen.TriplyPeriodicSurfaces import *
import numpy as np
import cadquery as cq

from OCP.StlAPI import StlAPI_Reader

from OCP.TopoDS import (
    TopoDS,
    TopoDS_Shape,
    TopoDS_Builder,
    TopoDS_Compound,
    TopoDS_Iterator,
    TopoDS_Wire,
    TopoDS_Face,
    TopoDS_Edge,
    TopoDS_Vertex,
    TopoDS_Solid,
    TopoDS_Shell,
    TopoDS_CompSolid,
)


#----------CYLINDER-----------------------------------------------------------------------------------------#

class tpms :
    def __init__(self,center,type_tpms,type_part,ske_type,n):
        self.center=center
        self.number=n
        self.name_part='tpms' + str(self.number)
        self.type_tpms=type_tpms
        self.type_part=type_part
        self.ske_type=ske_type

    def createSurfaces(self,path_data,rve,thickness,sizeMesh=0.05,minFacetAngle=10,maxRadius=0.05) :
        GenerateTPMS(self.type_tpms, thickness, rve, sizeMesh,minFacetAngle,maxRadius,path_data)

    def create_tpms(self,path_data,rve) :
    
        surf_tp = TopoDS_Shape()
        surf_tm = TopoDS_Shape()
        surf_p = TopoDS_Shape()
        surf_m = TopoDS_Shape()
        stl_reader = StlAPI_Reader()
        if path_data != '':
            stl_reader.Read(surf_tp, path_data + '/' + 'tpms_testplus.stl')
            stl_reader.Read(surf_tm, path_data + '/' + 'tpms_testminus.stl')
            stl_reader.Read(surf_p, path_data + '/' + 'tpms_plus.stl')
            stl_reader.Read(surf_m, path_data + '/' + 'tpms_minus.stl')
        else:
            stl_reader.Read(surf_tp, 'tpms_testplus.stl')
            stl_reader.Read(surf_tm, 'tpms_testminus.stl')
            stl_reader.Read(surf_p, 'tpms_plus.stl')
            stl_reader.Read(surf_m, 'tpms_minus.stl')

        face_cut_tp = cq.Face(surf_tp)
        face_cut_tm = cq.Face(surf_tm)
        face_cut_p = cq.Face(surf_p)
        face_cut_m = cq.Face(surf_m)

        box = cq.Workplane("front").box(rve.dx, rve.dy, rve.dz)
        
        boxCut = box.split(face_cut_p)
        boxCut = boxCut.split(face_cut_m)

        boxSolids = boxCut.solids().all()
        boxSolidsSize = boxCut.solids().size()
        
        print('boxSolids', boxSolids)
        print('boxSolidsSize', boxSolidsSize)
        print('boxCut', boxCut)

        listSolids = []

        for solid in boxSolids:
            temp = solid.split(face_cut_tp)
            temp = temp.split(face_cut_tm)
            tempSolids = temp.solids().all()
            tempSize = temp.solids().size()
            listSolids.append( (tempSize, solid.val()) )
#            print('tempSolids',tempSolids)
#            print('tempSize',tempSize)

        sheet = [el[1] for el in listSolids if el[0] > 1]
        skeletal = [el[1] for el in listSolids if el[0] == 1]

#        sheet = [el[1] for el in listSolids]
        
        print('sheet', sheet)
        print('skeletal', skeletal)

        if(self.type_part == 'sheet'):
            to_fuse = [cq.Shape(s.wrapped) for s in sheet]
            return_object = fuse_parts(to_fuse, True)
            return cq.Workplane().add(return_object[0])
        elif(self.type_part == 'skeletal'):
#            if(self.ske_type == 'plus'):
#                return cq.Workplane().add(final_pm)
#            elif(self.ske_type == 'minus'):
#                return cq.Workplane().add(final_mp)
#            elif(self.ske_type == 'double'):
            to_fuse = [cq.Shape(s.wrapped) for s in skeletal]
            return_object = fuse_parts(to_fuse, False)
            return cq.Workplane().add(return_object[0])
