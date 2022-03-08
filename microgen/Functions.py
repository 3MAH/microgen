import numpy as np
import cadquery as cq
import os

from OCP.BRepAlgoAPI import (
    BRepAlgoAPI_Common,
    BRepAlgoAPI_Fuse,
    BRepAlgoAPI_Cut,
    BRepAlgoAPI_BooleanOperation,
    BRepAlgoAPI_Splitter,
    BRepAlgoAPI_BuilderAlgo,
)

from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain


def rotateEuler(object,center,psi,theta,phi) :
    
    u = np.array([0., 0., 1.])
    u = np.array([np.cos(psi*np.pi/180.),np.sin(psi*np.pi/180.),0.])
    z2 = np.array([np.sin(psi*np.pi/180.)*np.sin(theta*np.pi/180.),-np.sin(theta*np.pi/180.)*np.cos(psi*np.pi/180.),np.cos(theta*np.pi/180.)])

    object_r = object.rotate((center[0],center[1],center[2]), (center[0],center[1],center[2]+1.0), psi)
    object_r = object_r.rotate((center[0],center[1],center[2]), (center[0]+u[0],center[1]+u[1],center[2]+u[2]), theta)
    object_r = object_r.rotate((center[0],center[1],center[2]), (center[0]+z2[0],center[1]+z2[1],center[2]+z2[2]), phi)
    return object_r

def remove_empty_lines(filename):
    if not os.path.isfile(filename):
        print("{} does not exist ".format(filename))
        return
    with open(filename) as filehandle:
        lines = filehandle.readlines()

    with open(filename, 'w') as filehandle:
        lines = filter(lambda x: x.strip(), lines)
        filehandle.writelines(lines)

def periodic(cqshape,rve):

    wk_plane = cq.Workplane().add(cqshape.Solids())
    periodic_object = []
    plist=[0,0,0,0,0,0]
    pcount=0
 
#    #création de vecteurs de translation des copies de test
#    tx=np.array([-rve.dx/2.0,3.0*rve.dx/2.0,rve.dx/2.0,rve.dx/2.0,rve.dx/2.0,rve.dx/2.0])
#    ty=np.array([rve.dy/2.0,rve.dy/2.0,-rve.dy/2.0,3.0*rve.dy/2.0,rve.dy/2.0,rve.dy/2.0])
#    tz=np.array([rve.dz/2.0,rve.dz/2.0,rve.dz/2.0,rve.dz/2.0,-rve.dz/2.0,3.0*rve.dz/2.0])
#
#    #création des vecteurs de translation des copies définitives
#    x_faces = {'1':0.,'2':rve.dx,'3':0.,'4':0.,'5':0.,'6':0.}
#    y_faces = {'1':0.,'2':0.,'3':0.,'4':rve.dy,'5':0.,'6':0.}
#    z_faces = {'1':0.,'2':0.,'3':0,'4':0.,'5':0.,'6':rve.dz}
#
#    x_edges = {'13':0.,'14':0.,'15':0,'16':0,'23':rve.dx,'24':rve.dx,'25':rve.dx,'26':rve.dx,'35':0.,'36':0.,'45':0.,'46':0.}
#    y_edges = {'13':0.,'14':rve.dy,'15':0.,'16':0.,'23':0,'24':rve.dy,'25':0.,'26':0.,'35':0,'36':0.,'45':rve.dy,'46':rve.dy}
#    z_edges = {'13':0.,'14':0.,'15':0.,'16':rve.dz,'23':0.,'24': 0.,'25':0.,'26':rve.dz,'35':0,'36':rve.dz,'45':0.,'46': rve.dz}
#
#    x_corners = {'135':0.,'136':0.,'145':0,'146':0,'235':rve.dx,'236':rve.dx,'245':rve.dx,'246':rve.dx}
#    y_corners = {'135':0.,'136':0.,'145':rve.dy,'146':rve.dz,'235':0.,'236':0.,'245':rve.dy,'246':rve.dx}
#    z_corners = {'135':0.,'136':rve.dz,'145':0,'146':rve.dz,'235':0.,'236':rve.dz,'245':0.,'246':rve.dx}

    #creation de faces de coupe
    Plane_xm = cq.Face.makePlane(basePnt = (0, 0, 0), dir = (-1, 0, 0))
    Plane_xp = cq.Face.makePlane(basePnt = (rve.dx, 0, 0), dir = (1, 0, 0))
    Plane_ym = cq.Face.makePlane(basePnt = (0, 0, 0), dir = (0, -1, 0))
    Plane_yp = cq.Face.makePlane(basePnt = (0, rve.dy, 0), dir = (0, 1, 0))
    Plane_zm = cq.Face.makePlane(basePnt = (0, 0, 0), dir = (0, 0, -1))
    Plane_zp = cq.Face.makePlane(basePnt = (0, 0, rve.dz), dir = (0, 0, 1))
#
    #Test des partitions
    Partition_xm = wk_plane.split(cq.Workplane().add(Plane_xm))
    Partition_xp = wk_plane.split(cq.Workplane().add(Plane_xp))
    Partition_ym = wk_plane.split(cq.Workplane().add(Plane_ym))
    Partition_yp = wk_plane.split(cq.Workplane().add(Plane_yp))
    Partition_zm = wk_plane.split(cq.Workplane().add(Plane_zm))
    Partition_zp = wk_plane.split(cq.Workplane().add(Plane_zp))
    
    Partition=[Partition_xm.solids().all(), Partition_xp.solids().all(), Partition_ym.solids().all(), Partition_yp.solids().all(), Partition_zm.solids().all(), Partition_zp.solids().all()]

    for i in range (0,6):
        if len(Partition[i]) > 1:
            plist[i]=1
            pcount=pcount+1

    for i in range(0,3):
        if plist[i*2] + plist[i*2+1] == 2:
            plist[i*2] = 0
            pcount=pcount-1
            plist[i*2+1] = 0
            pcount=pcount-1

    print(plist)

    if pcount==0: #inclusion ne dépasse pas du cube
        periodic_object.append(cqshape.solids().intersect(rve.Box))
    elif pcount==1:
        if plist[0] == 1:
            periodic_object.append(Partition_xm.solids(">X").intersect(rve.Box))
            periodic_object.append(Partition_xm.solids("<X").translate((rve.dx,0,0)).intersect(rve.Box))
        if plist[1] == 1:
            periodic_object.append(Partition_xp.solids("<X").intersect(rve.Box))
            periodic_object.append(Partition_xp.solids(">X").translate((-rve.dx,0,0)).intersect(rve.Box))
        if plist[2] == 1:
            periodic_object.append(Partition_ym.solids(">Y").intersect(rve.Box))
            periodic_object.append(Partition_ym.solids("<Y").translate((0,rve.dy,0)).intersect(rve.Box))
        if plist[3] == 1:
            periodic_object.append(Partition_yp.solids("<Y").intersect(rve.Box))
            periodic_object.append(Partition_yp.solids(">Y").translate((0,-rve.dy,0)).intersect(rve.Box))
        if plist[4] == 1:
            periodic_object.append(Partition_zm.solids(">Z").intersect(rve.Box))
            periodic_object.append(Partition_zm.solids("<Z").translate((0,0,rve.dz)).intersect(rve.Box))
        if plist[5] == 1:
            periodic_object.append(Partition_zp.solids("<Z").intersect(rve.Box))
            periodic_object.append(Partition_zp.solids(">Z").translate((0,0,-rve.dz)).intersect(rve.Box))
    elif pcount==2:
        ed=''
        for k in range(len(plist)):
            if plist[k]==1:
                ed+=str(k+1)
        edge=int(ed)
        
        if edge==13:
            Partition_xpy = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym))
            periodic_object.append(Partition_xpy.solids(">Y").intersect(rve.Box))
            periodic_object.append(Partition_xpy.solids("<Y").translate((0,rve.dy,0)).intersect(rve.Box))
            Partition_xmy = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym))
            periodic_object.append(Partition_xmy.solids(">Y").translate((rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmy.solids("<Y").translate((rve.dx,rve.dy,0)).intersect(rve.Box))
        if edge==14:
            Partition_xpy = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_yp))
            periodic_object.append(Partition_xpy.solids("<Y").intersect(rve.Box))
            periodic_object.append(Partition_xpy.solids(">Y").translate((0,-rve.dy,0)).intersect(rve.Box))
            Partition_xmy = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_yp))
            periodic_object.append(Partition_xmy.solids("<Y").translate((rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmy.solids(">Y").translate((rve.dx,-rve.dy,0)).intersect(rve.Box))
        if edge==23:
            Partition_xmy = cq.Workplane().add(Partition_xp.solids("<X")).split(cq.Workplane().add(Plane_ym))
            periodic_object.append(Partition_xmy.solids(">Y").intersect(rve.Box))
            periodic_object.append(Partition_xmy.solids("<Y").translate((0,rve.dy,0)).intersect(rve.Box))
            Partition_xpy = cq.Workplane().add(Partition_xp.solids(">X")).split(cq.Workplane().add(Plane_ym))
            periodic_object.append(Partition_xpy.solids(">Y").translate((-rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xpy.solids("<Y").translate((-rve.dx,rve.dy,0)).intersect(rve.Box))
        if edge==24:
            Partition_xmy = cq.Workplane().add(Partition_xp.solids("<X")).split(cq.Workplane().add(Plane_yp))
            periodic_object.append(Partition_xmy.solids("<Y").intersect(rve.Box))
            periodic_object.append(Partition_xmy.solids(">Y").translate((0,-rve.dy,0)).intersect(rve.Box))
            Partition_xpy = cq.Workplane().add(Partition_xp.solids(">X")).split(cq.Workplane().add(Plane_yp))
            periodic_object.append(Partition_xpy.solids("<Y").translate((-rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xpy.solids(">Y").translate((-rve.dx,-rve.dy,0)).intersect(rve.Box))
            
        if edge==15:
            Partition_xpz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_zm))
            periodic_object.append(Partition_xpz.solids(">Z").intersect(rve.Box))
            periodic_object.append(Partition_xpz.solids("<Z").translate((0,0,rve.dz)).intersect(rve.Box))
            Partition_xmz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_zm))
            periodic_object.append(Partition_xmz.solids(">Z").translate((rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmz.solids("<Z").translate((rve.dx,0,rve.dz)).intersect(rve.Box))
        if edge==16:
            Partition_xpz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_zp))
            periodic_object.append(Partition_xpz.solids("<Z").intersect(rve.Box))
            periodic_object.append(Partition_xpz.solids(">Z").translate((0,0,-rve.dz)).intersect(rve.Box))
            Partition_xmz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_zp))
            periodic_object.append(Partition_xmz.solids("<Z").translate((rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmz.solids(">Z").translate((rve.dx,0,-rve.dz)).intersect(rve.Box))
        if edge==25:
            Partition_xmz = cq.Workplane().add(Partition_xp.solids("<X")).split(cq.Workplane().add(Plane_zm))
            periodic_object.append(Partition_xmz.solids(">Z").intersect(rve.Box))
            periodic_object.append(Partition_xmz.solids("<Z").translate((0,0,rve.dz)).intersect(rve.Box))
            Partition_xpz = cq.Workplane().add(Partition_xp.solids(">X")).split(cq.Workplane().add(Plane_zm))
            periodic_object.append(Partition_xpz.solids(">Z").translate((-rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xpz.solids("<Z").translate((-rve.dx,0,rve.dz)).intersect(rve.Box))
        if edge==26:
            Partition_xmz = cq.Workplane().add(Partition_xp.solids("<X")).split(cq.Workplane().add(Plane_zp))
            periodic_object.append(Partition_xmz.solids("<Z").intersect(rve.Box))
            periodic_object.append(Partition_xmz.solids(">Z").translate((0,0,-rve.dz)).intersect(rve.Box))
            Partition_xpz = cq.Workplane().add(Partition_xp.solids(">X")).split(cq.Workplane().add(Plane_zp))
            periodic_object.append(Partition_xpz.solids("<Z").translate((-rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xpz.solids(">Z").translate((-rve.dx,0,-rve.dz)).intersect(rve.Box))
            
        if edge==35:
            Partition_ypz = cq.Workplane().add(Partition_ym.solids(">Y")).split(cq.Workplane().add(Plane_zm))
            periodic_object.append(Partition_ypz.solids(">Z").intersect(rve.Box))
            periodic_object.append(Partition_ypz.solids("<Z").translate((0,0,rve.dz)).intersect(rve.Box))
            Partition_ymz = cq.Workplane().add(Partition_ym.solids("<Y")).split(cq.Workplane().add(Plane_zm))
            periodic_object.append(Partition_ymz.solids(">Z").translate((0,rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_ymz.solids("<Z").translate((0,rve.dy,rve.dz)).intersect(rve.Box))
        if edge==36:
            Partition_ypz = cq.Workplane().add(Partition_ym.solids(">Y")).split(cq.Workplane().add(Plane_zp))
            periodic_object.append(Partition_ypz.solids("<Z").intersect(rve.Box))
            periodic_object.append(Partition_ypz.solids(">Z").translate((0,0,-rve.dz)).intersect(rve.Box))
            Partition_ymz = cq.Workplane().add(Partition_ym.solids("<Y")).split(cq.Workplane().add(Plane_zp))
            periodic_object.append(Partition_ymz.solids("<Z").translate((0,rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_ymz.solids(">Z").translate((0,rve.dy,-rve.dz)).intersect(rve.Box))
        if edge==45:
            Partition_ymz = cq.Workplane().add(Partition_yp.solids("<Y")).split(cq.Workplane().add(Plane_zm))
            periodic_object.append(Partition_ymz.solids(">Z").intersect(rve.Box))
            periodic_object.append(Partition_ymz.solids("<Z").translate((0,0,rve.dz)).intersect(rve.Box))
            Partition_ypz = cq.Workplane().add(Partition_yp.solids(">Y")).split(cq.Workplane().add(Plane_zm))
            periodic_object.append(Partition_ypz.solids(">Z").translate((0,-rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_ypz.solids("<Z").translate((0,-rve.dy,rve.dz)).intersect(rve.Box))
        if edge==46:
            Partition_ymz = cq.Workplane().add(Partition_yp.solids("<Y")).split(cq.Workplane().add(Plane_zp))
            periodic_object.append(Partition_ymz.solids("<Z").intersect(rve.Box))
            periodic_object.append(Partition_ymz.solids(">Z").translate((0,0,-rve.dz)).intersect(rve.Box))
            Partition_ypz = cq.Workplane().add(Partition_yp.solids(">Y")).split(cq.Workplane().add(Plane_zp))
            periodic_object.append(Partition_ypz.solids("<Z").translate((0,-rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_ypz.solids(">Z").translate((0,-rve.dy,-rve.dz)).intersect(rve.Box))

    elif pcount==3: #inclusion sur un coin
        #détection du coin concerné
        co=''
        for k in range(len(plist)):
            if plist[k]==1:
                co+=str(k+1)
        corner=int(co)

        if corner==135:
            Partition_xpypz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xpypz.solids(">Z").intersect(rve.Box))
            periodic_object.append(Partition_xpypz.solids("<Z").translate((0,0,rve.dz)).intersect(rve.Box))
            Partition_xpymz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xpymz.solids(">Z").translate((0,rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xpymz.solids("<Z").translate((0,rve.dy,rve.dz)).intersect(rve.Box))
            Partition_xmypz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xmypz.solids(">Z").translate((rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmypz.solids("<Z").translate((rve.dx,0,rve.dz)).intersect(rve.Box))
            Partition_xmymz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xmymz.solids(">Z").translate((rve.dx,rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmymz.solids("<Z").translate((rve.dx,rve.dy,rve.dz)).intersect(rve.Box))
        if corner==136:
            Partition_xpypz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xpypz.solids("<Z").intersect(rve.Box))
            periodic_object.append(Partition_xpypz.solids(">Z").translate((0,0,-rve.dz)).intersect(rve.Box))
            Partition_xpymz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xpymz.solids("<Z").translate((0,rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xpymz.solids(">Z").translate((0,rve.dy,-rve.dz)).intersect(rve.Box))
            Partition_xmypz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xmypz.solids("<Z").translate((rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmypz.solids(">Z").translate((rve.dx,0,-rve.dz)).intersect(rve.Box))
            Partition_xmymz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xmymz.solids("<Z").translate((rve.dx,rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmymz.solids(">Z").translate((rve.dx,rve.dy,-rve.dz)).intersect(rve.Box))
        if corner==145:
            Partition_xpypz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xpypz.solids(">Z").intersect(rve.Box))
            periodic_object.append(Partition_xpypz.solids("<Z").translate((0,0,rve.dz)).intersect(rve.Box))
            Partition_xpymz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xpymz.solids(">Z").translate((0,-rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xpymz.solids("<Z").translate((0,-rve.dy,rve.dz)).intersect(rve.Box))
            Partition_xmypz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xmypz.solids(">Z").translate((rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmypz.solids("<Z").translate((rve.dx,0,rve.dz)).intersect(rve.Box))
            Partition_xmymz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xmymz.solids(">Z").translate((rve.dx,-rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmymz.solids("<Z").translate((rve.dx,-rve.dy,rve.dz)).intersect(rve.Box))
        if corner==146:
            Partition_xpypz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xpypz.solids("<Z").intersect(rve.Box))
            periodic_object.append(Partition_xpypz.solids(">Z").translate((0,0,-rve.dz)).intersect(rve.Box))
            Partition_xpymz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xpymz.solids("<Z").translate((0,-rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xpymz.solids(">Z").translate((0,-rve.dy,-rve.dz)).intersect(rve.Box))
            Partition_xmypz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xmypz.solids("<Z").translate((rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmypz.solids(">Z").translate((rve.dx,0,-rve.dz)).intersect(rve.Box))
            Partition_xmymz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xmymz.solids("<Z").translate((rve.dx,-rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmymz.solids(">Z").translate((rve.dx,-rve.dy,-rve.dz)).intersect(rve.Box))
        if corner==235:
            Partition_xpypz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xpypz.solids(">Z").intersect(rve.Box))
            periodic_object.append(Partition_xpypz.solids("<Z").translate((0,0,rve.dz)).intersect(rve.Box))
            Partition_xpymz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xpymz.solids(">Z").translate((0,rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xpymz.solids("<Z").translate((0,rve.dy,rve.dz)).intersect(rve.Box))
            Partition_xmypz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xmypz.solids(">Z").translate((-rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmypz.solids("<Z").translate((-rve.dx,0,rve.dz)).intersect(rve.Box))
            Partition_xmymz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xmymz.solids(">Z").translate((-rve.dx,rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmymz.solids("<Z").translate((-rve.dx,rve.dy,rve.dz)).intersect(rve.Box))
        if corner==236:
            Partition_xpypz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xpypz.solids("<Z").intersect(rve.Box))
            periodic_object.append(Partition_xpypz.solids(">Z").translate((0,0,-rve.dz)).intersect(rve.Box))
            Partition_xpymz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xpymz.solids("<Z").translate((0,rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xpymz.solids(">Z").translate((0,rve.dy,-rve.dz)).intersect(rve.Box))
            Partition_xmypz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xmypz.solids("<Z").translate((-rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmypz.solids(">Z").translate((-rve.dx,0,-rve.dz)).intersect(rve.Box))
            Partition_xmymz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xmymz.solids("<Z").translate((-rve.dx,rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmymz.solids(">Z").translate((-rve.dx,rve.dy,-rve.dz)).intersect(rve.Box))
        if corner==245:
            Partition_xpypz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xpypz.solids(">Z").intersect(rve.Box))
            periodic_object.append(Partition_xpypz.solids("<Z").translate((0,0,rve.dz)).intersect(rve.Box))
            Partition_xpymz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xpymz.solids(">Z").translate((0,-rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xpymz.solids("<Z").translate((0,-rve.dy,rve.dz)).intersect(rve.Box))
            Partition_xmypz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xmypz.solids(">Z").translate((-rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmypz.solids("<Z").translate((-rve.dx,0,rve.dz)).intersect(rve.Box))
            Partition_xmymz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xmymz.solids(">Z").translate((-rve.dx,-rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmymz.solids("<Z").translate((-rve.dx,-rve.dy,rve.dz)).intersect(rve.Box))
        if corner==246:
            Partition_xpypz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xpypz.solids("<Z").intersect(rve.Box))
            periodic_object.append(Partition_xpypz.solids(">Z").translate((0,0,-rve.dz)).intersect(rve.Box))
            Partition_xpymz = cq.Workplane().add(Partition_xm.solids("<X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xpymz.solids("<Z").translate((0,-rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xpymz.solids(">Z").translate((0,-rve.dy,-rve.dz)).intersect(rve.Box))
            Partition_xmypz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids("<Y")
            periodic_object.append(Partition_xmypz.solids("<Z").translate((-rve.dx,0,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmypz.solids(">Z").translate((-rve.dx,0,-rve.dz)).intersect(rve.Box))
            Partition_xmymz = cq.Workplane().add(Partition_xm.solids(">X")).split(cq.Workplane().add(Plane_ym)).solids(">Y")
            periodic_object.append(Partition_xmymz.solids("<Z").translate((-rve.dx,-rve.dy,0)).intersect(rve.Box))
            periodic_object.append(Partition_xmymz.solids(">Z").translate((-rve.dx,-rve.dy,-rve.dz)).intersect(rve.Box))

    
    occ_solids_list = [s.val().Solids() for s in periodic_object]
    flat_list = [item.copy() for sublist in occ_solids_list for item in sublist]
    to_fuse = [cq.Shape(s.wrapped) for s in flat_list]
    return_object_periodic = fuse_parts(to_fuse, False)
    return (return_object_periodic[0].copy(), flat_list)

def fuse_parts(CqShapeList, retain_edges):
          
#    occ_solids_list = (s.Solids() for s in CqShapeList)
#    for cqshape in CqShapeList:

    occ_solids_list = [s.Solids() for s in CqShapeList]
#    print("occ_solids_list = ", occ_solids_list)
#    flat_list = [item for sublist in occ_solids_list for item in sublist]
    
 #   print("flat_list = ", flat_list)

    occ_Solids = CqShapeList[0].wrapped
    for i in range(1, len(CqShapeList)):
        fuse = BRepAlgoAPI_Fuse(occ_Solids, CqShapeList[i].wrapped)
        occ_Solids = fuse.Shape()

    if retain_edges == True:
        return (cq.Shape(occ_Solids), occ_solids_list)
    else:
        upgrader = ShapeUpgrade_UnifySameDomain(occ_Solids, True, True, True)
        upgrader.Build()
        fixed = upgrader.Shape()
        occ_solids_list = [[cq.Solid(fixed)]]

        return (cq.Shape(fixed), occ_solids_list)

#def cut_parts(CqShapeList):
#
#    print('inside cut')
#    phase_cut = []
#    phase_cut.append(CqShapeList[0].copy())
#    cut_objtemp = CqShapeList[0].copy()
#    upgrader = ShapeUpgrade_UnifySameDomain(cut_objtemp.wrapped, True, True, True)
#    upgrader.Build()
#    cut_obj = cq.Shape(upgrader.Shape())
#
#    for shape in CqShapeList[1::]:
#        print('tatayoyo')
#
#        SolidsCut = []
#        for s in shape.Solids():
#            sCut = s.wrapped
#            for t in cut_obj.Solids():
#                cut = BRepAlgoAPI_Cut(sCut, t.wrapped)
#                sCut = cut.Shape()
#            SolidsCut.append(cq.Shape(cut.Shape()))
#        cutted = fuse_parts(SolidsCut, False)
#        phase_cut.append(cutted[0])
#
#        fuse = BRepAlgoAPI_Fuse(cut_obj.wrapped, shape.wrapped)
#        fused = fuse.Shape()
#        upgrader = ShapeUpgrade_UnifySameDomain(fused, True, True, True)
#        upgrader.Build()
#        cut_obj = cq.Shape(upgrader.Shape())
#
#    occ_solids_list = [s.Solids() for s in phase_cut]
#    print(phase_cut)
#    print(occ_solids_list)
#    print('outside cut')
#
#    return (phase_cut, occ_solids_list)

def cutPhasesByShape(CqShapeList, cut_obj):
    phase_cut = []

    for shape in CqShapeList:
        cut = BRepAlgoAPI_Cut(shape.wrapped, cut_obj.wrapped)
        if len(cq.Shape(cut.Shape()).Solids()) > 0:
            phase_cut.append(cq.Shape(cut.Shape()))

    occ_solids_list = [s.Solids() for s in phase_cut]
    print(phase_cut)
    print(occ_solids_list)
    print('outside cut')

    return (phase_cut, occ_solids_list)

def cutPhaseByShapeList(phaseToCut, CqShapeList):

    ResultCut = phaseToCut
    for shape in CqShapeList:
        cut = BRepAlgoAPI_Cut(ResultCut.wrapped, shape.wrapped)
        ResultCut = cq.Shape(cut.Shape())

    occ_solids_list = ResultCut.Solids()
    return (ResultCut, occ_solids_list)


def cut_parts(CqShapeList, reverseOrder = True):
    phase_cut = []
    if(reverseOrder == True):
        CqShapeList_inv = CqShapeList[::-1]
    else:
        CqShapeList_inv = CqShapeList
#    print(CqShapeList)
#    print(CqShapeList_inv)

    cut_obj = CqShapeList_inv[0].copy()
    phase_cut.append(cut_obj)

    for shape in CqShapeList_inv[1::]:
        copy = shape.copy()
        cut = BRepAlgoAPI_Cut(copy.wrapped, cut_obj.wrapped)
        phase_cut.append(cq.Shape(cut.Shape()))

        fuse = BRepAlgoAPI_Fuse(cut_obj.wrapped, shape.wrapped)
        fused = fuse.Shape()
        upgrader = ShapeUpgrade_UnifySameDomain(fused, True, True, True)
        upgrader.Build()
        cut_obj = cq.Shape(upgrader.Shape())

    phase_cut.reverse()

    occ_solids_list = [s.Solids() for s in phase_cut]
    print(phase_cut)
    print(occ_solids_list)
    print('outside cut')

    return (phase_cut, occ_solids_list)


#def cut_parts(CqShapeList):
#
#    phase_cut = []
#    occ_Solids = CqShapeList[-1].copy()
#    phase_cut.append(CqShapeList[-1])
#
#    for s in CqShapeList[-2::-1]:
#        print('s', s)
#        cut = BRepAlgoAPI_Cut(s.wrapped, occ_Solids.wrapped)
#        phase_cut.append(cq.Shape(cut.Shape()))
#
#        fuse = BRepAlgoAPI_Fuse(occ_Solids.wrapped, s.wrapped)
#        occ_Solids = fuse.Shape()
#        upgrader = ShapeUpgrade_UnifySameDomain(occ_Solids, True, True, True)
#        upgrader.Build()
#        occ_Solids = cq.Shape(upgrader.Shape())
#
#    print('phase_cut', phase_cut)
#    occ_solids_list = [s.Solids() for s in phase_cut[::-1]]
#    return (phase_cut[::-1], occ_solids_list)
    
# Ajout MB 07/01/2022

def lanceNeper(filename,nbCell,dimCube):
    command = "neper -T -n "+str(nbCell)+ " -id 1 -dim 3"
    command = command + " -domain 'cube(" + str(dimCube[0]) + "," + str(dimCube[1]) + "," + str(dimCube[2]) + ")'"
    command = command + " -morpho ""gg"" -o " + filename
    
    os.system(command)
    

def parseNeper(filename):
    
    # 1ere étape : lecture des coordonées des germes, coordonées des sommets,
    # labels globaux des sommets pour chaque segment, labels globaux des segments
    # pour chaque face, labels globaux des faces pour chaque polyèdre à partir
    # du fichier de tesselation filename.tess généré par NEPER
    
    file=filename + '.tess'
    fid = open(file,'r')
    
    flagChercheSeed = False
    seed = []
    
    flagChercheVertex = False
    vertices = []
    
    flagChercheEdges = False
    edges = []
    
    flagChercheFaces = False
    faces = []
    iLigneFace = 1
    
    flagCherchePolys = False
    polys = []
    
    flagExtraction = False
    i=0
    

    for line in fid:
        
        if (flagExtraction):
            if '*' in line:
                flagExtraction = False
            else:
                if i == 1: # Extraction des germes
                    ligneDecoupee = line.split()
                    seed.append([float(ligneDecoupee[1]),float(ligneDecoupee[2]),float(ligneDecoupee[3])])
                if i == 2: # Extraction des sommets
                    ligneDecoupee = line.split()
                    vertices.append([int(ligneDecoupee[0]),float(ligneDecoupee[1]),float(ligneDecoupee[2]),float(ligneDecoupee[3])])
                if i == 3: # Extraction des segments
                    ligneDecoupee = line.split()
                    edges.append([int(ligneDecoupee[0]),int(ligneDecoupee[1]),int(ligneDecoupee[2])])
                if i == 4: # Extraction des faces
                    if iLigneFace == 1:
                        tmp = []
                        tmp.append(int(line.split()[0]))
                    if iLigneFace == 2:
                        ligneDecoupee = line.split()
                        nbEdges = int(ligneDecoupee[0])
                        for k in range(nbEdges):
                            tmp.append(int(ligneDecoupee[k+1]))
                        faces.append(tmp)
                    if iLigneFace == 4:
                        iLigneFace = 0
                        
                    iLigneFace +=1
                if i == 5: # Extraction des cellules
                    tmp = []
                    ligneDecoupee = line.split()
                    tmp.append(int(ligneDecoupee[0]))
                    nbFaces = int(ligneDecoupee[1])
                    for k in range(nbFaces):
                        tmp.append(np.abs(int(ligneDecoupee[k+2])))
                    polys.append(tmp)
                
                        

        if '*seed' in line:
            flagChercheSeed = True
            i+=1
            
        if '**vertex' in line:
            flagChercheVertex = True
            i+=1
            next(fid)

        if '**edge' in line:
            flagChercheEdges = True
            i+=1
            next(fid)
            
        if '**face' in line:
            flagChercheFaces = True
            i+=1
            next(fid)
            
        if '**polyhedron' in line:
            flagCherchePolys = True
            i+=1
            next(fid)
            
        if (flagChercheSeed):
            flagExtraction = True
            flagChercheSeed = False
            
        if (flagChercheVertex):
            flagExtraction = True
            flagChercheVertex = False
            
        if (flagChercheEdges):
            flagExtraction = True
            flagChercheEdges = False

        if (flagChercheFaces):
            flagExtraction = True
            flagChercheFaces = False

        if (flagCherchePolys):
            flagExtraction = True
            flagCherchePolys = False
            
    fid.close()

    # 2e étape : création d'un objet pour le polycristal ayant une structure semblable
    # à la sortie de compute_voronoi de la librairie py_voro, compatible avec la définition
    # des polyèdres de la librairie Microgen

    A = []
    nbPolys = len(polys)
    listeSommets=[]
    listeSommetsOut = []
    for i in range(nbPolys):
        voro = {}
        voro["original"] = seed[i]
        voro["faces"]=[]
        listeSommetsPoly = []
        sommets=[]
        for facePoly in polys[i][1:]:
            listeSommetsFace=[]
            dicVertices = {}
            dicVertices["vertices"]=[]
            for segment in faces[facePoly-1][1:]:
                # Listes des sommets avec numérotation globale et faces associées
                if edges[np.abs(segment)-1][1] not in listeSommets:
                    dicFacesAssociees = {}
                    dicFacesAssociees["cell_associees"]=[facePoly]
                    listeSommets.append(edges[np.abs(segment)-1][1])
                    listeSommetsOut.append((edges[np.abs(segment)-1][1],dicFacesAssociees))
                else:
                    idx = listeSommets.index(edges[np.abs(segment)-1][1])
                    if facePoly not in listeSommetsOut[idx][1]["cell_associees"]:
                        listeSommetsOut[idx][1]["cell_associees"].append(facePoly)
                    #sommets.append(vertices[edges[np.abs(segment)-1][1]-1][1:])
                if edges[np.abs(segment)-1][2] not in listeSommets:
                    dicFacesAssociees = {}
                    dicFacesAssociees["cell_associees"]=[facePoly]
                    listeSommets.append(edges[np.abs(segment)-1][2])
                    listeSommetsOut.append((edges[np.abs(segment)-1][2],dicFacesAssociees))
                else:
                    idx = listeSommets.index(edges[np.abs(segment)-1][2])
                    if facePoly not in listeSommetsOut[idx][1]["cell_associees"]:
                        listeSommetsOut[idx][1]["cell_associees"].append(facePoly)
                    #sommets.append(vertices[edges[np.abs(segment)-1][2]-1][1:])
                
                # Listes des sommets et de leurs coordonnées avec numérotation locale
                # à inclure dans voro
                if segment < 0:
                    if edges[np.abs(segment)-1][2] not in listeSommetsPoly:
                        listeSommetsPoly.append(edges[np.abs(segment)-1][2])
                        sommets.append(vertices[edges[np.abs(segment)-1][2]-1][1:])
                    if edges[np.abs(segment)-1][2] not in listeSommetsFace:
                        listeSommetsFace.append(edges[np.abs(segment)-1][2])
                        dicVertices["vertices"].append(listeSommetsPoly.index(edges[np.abs(segment)-1][2]))
                    
                    if edges[np.abs(segment)-1][1] not in listeSommetsPoly:
                        listeSommetsPoly.append(edges[np.abs(segment)-1][1])
                        sommets.append(vertices[edges[np.abs(segment)-1][1]-1][1:])
                    if edges[np.abs(segment)-1][1] not in listeSommetsFace:
                        listeSommetsFace.append(edges[np.abs(segment)-1][1])
                        dicVertices["vertices"].append(listeSommetsPoly.index(edges[np.abs(segment)-1][1]))

                if segment > 0:
                    if edges[np.abs(segment)-1][1] not in listeSommetsPoly:
                        listeSommetsPoly.append(edges[np.abs(segment)-1][1])
                        sommets.append(vertices[edges[np.abs(segment)-1][1]-1][1:])
                    if edges[np.abs(segment)-1][1] not in listeSommetsFace:
                        listeSommetsFace.append(edges[np.abs(segment)-1][1])
                        dicVertices["vertices"].append(listeSommetsPoly.index(edges[np.abs(segment)-1][1]))
                    
                    if edges[np.abs(segment)-1][2] not in listeSommetsPoly:
                        listeSommetsPoly.append(edges[np.abs(segment)-1][2])
                        sommets.append(vertices[edges[np.abs(segment)-1][2]-1][1:])
                    if edges[np.abs(segment)-1][2] not in listeSommetsFace:
                        listeSommetsFace.append(edges[np.abs(segment)-1][2])
                        dicVertices["vertices"].append(listeSommetsPoly.index(edges[np.abs(segment)-1][2]))
            voro["faces"].append(dicVertices)
                    
        voro["vertices"]=sommets
        A.append(voro)
        

             
    return A, seed, listeSommetsOut, edges, faces, polys

# fin ajout
