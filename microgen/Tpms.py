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
    
        surf_c = TopoDS_Shape()
        surf_p = TopoDS_Shape()
        surf_m = TopoDS_Shape()
        stl_reader = StlAPI_Reader()
        if path_data != '':
#            stl_reader.Read(surf_c, path_data + '/' + 'tpms_center.stl')
            stl_reader.Read(surf_p, path_data + '/' + 'tpms_plus.stl')
            stl_reader.Read(surf_m, path_data + '/' + 'tpms_minus.stl')
        else:
 #           stl_reader.Read(surf_c, 'tpms_center.stl')
            stl_reader.Read(surf_p, 'tpms_plus.stl')
            stl_reader.Read(surf_m, 'tpms_minus.stl')

#        face_cut = cq.Face(test)
        face_cut_p = cq.Face(surf_p)
        face_cut_m = cq.Face(surf_m)

        box = cq.Workplane("front").box(rve.dx, rve.dy, rve.dz)
        #box = box.translate((0.75,0.75,0.75))
        final_p = box.split(face_cut_p).solids(">X")
        final_pm = final_p.split(face_cut_m).solids("<X").val()
        final_pp = final_p.split(face_cut_m).solids(">X").val()

        final_m = box.split(face_cut_m).solids("<X")
        final_mp = final_m.split(face_cut_p).solids(">X").val()
        final_mm = final_m.split(face_cut_p).solids("<X").val()

        final_c = box.split(face_cut_p).solids("<X")
        final_c = final_c.split(face_cut_m).solids(">X").val()

        if(self.type_part == 'sheet'):
#            to_fuse = [cq.Shape(final_c.wrapped), cq.Shape(final_pp.wrapped), cq.Shape(final_mm.wrapped)]
#            return_object = fuse_parts(to_fuse, False)
#            return cq.Workplane().add(return_object[0])
            return cq.Workplane().add(final_c)
        elif(self.type_part == 'skeletal'):
            if(self.ske_type == 'plus'):
                return cq.Workplane().add(final_pm)
            elif(self.ske_type == 'minus'):
                return cq.Workplane().add(final_mp)
            elif(self.ske_type == 'double'):
                to_fuse = [cq.Shape(final_pm.wrapped), cq.Shape(final_mp.wrapped)]
                return_object = fuse_parts(to_fuse, False)
                return cq.Workplane().add(return_object[0])
