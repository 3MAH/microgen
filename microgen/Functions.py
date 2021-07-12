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

def periodic(object,rve):
	
    periodic_object = []
 
    #création de vecteurs de translation des copies de test
    tx=np.array([-rve.dx/2.0,3.0*rve.dx/2.0,rve.dx/2.0,rve.dx/2.0,rve.dx/2.0,rve.dx/2.0])
    ty=np.array([rve.dy/2.0,rve.dy/2.0,-rve.dy/2.0,3.0*rve.dy/2.0,rve.dy/2.0,rve.dy/2.0])
    tz=np.array([rve.dz/2.0,rve.dz/2.0,rve.dz/2.0,rve.dz/2.0,-rve.dz/2.0,3.0*rve.dz/2.0])
	
	#parcours des inclusions
    if object.shape!='matrix':
        plist=[0,0,0,0,0,0]
        phase_perios_test = []

	#création des vecteurs de translation des copies définitives
    x_faces = {'1':0.,'2':rve.dx,'3':0.,'4':0.,'5':0.,'6':0.}
    y_faces = {'1':0.,'2':0.,'3':0.,'4':rve.dy,'5':0.,'6':0.}
    z_faces = {'1':0.,'2':0.,'3':0,'4':0.,'5':0.,'6':rve.dz}
	
    x_edges = {'13':0.,'14':0.,'15':0,'16':0,'23':rve.dx,'24':rve.dx,'25':rve.dx,'26':rve.dx,'35':0.,'36':0.,'45':0.,'46':0.}
    y_edges = {'13':0.,'14':rve.dy,'15':0.,'16':0.,'23':0,'24':rve.dy,'25':0.,'26':0.,'35':0,'36':0.,'45':rve.dy,'46':rve.dy}
    z_edges = {'13':0.,'14':0.,'15':0.,'16':rve.dz,'23':0.,'24': 0.,'25':0.,'26':rve.dz,'35':0,'36':rve.dz,'45':0.,'46': rve.dz}
		
    x_corners = {'135':0.,'136':0.,'145':0,'146':0,'235':rve.dx,'236':rve.dx,'245':rve.dx,'246':rve.dx}
    y_corners = {'135':0.,'136':0.,'145':rve.dy,'146':rve.dz,'235':0.,'236':0.,'245':rve.dy,'246':rve.dx}
    z_corners = {'135':0.,'136':rve.dz,'145':0,'146':rve.dz,'235':0.,'236':rve.dz,'245':0.,'246':rve.dx}

	#creation de faces de coupe
    Plane_xm = cq.Face.makePlane(basePnt = (0, 0, 0), dir = (-1, 0, 0))
    Plane_xp = cq.Face.makePlane(basePnt = (rve.dx, 0, 0), dir = (1, 0, 0))
    Plane_ym = cq.Face.makePlane(basePnt = (0, 0, 0), dir = (0, -1, 0))
    Plane_yp = cq.Face.makePlane(basePnt = (0, rve.dy, 0), dir = (0, 1, 0))
    Plane_zm = cq.Face.makePlane(basePnt = (0, 0, 0), dir = (0, 0, -1))
    Plane_zp = cq.Face.makePlane(basePnt = (0, 0, rve.dz), dir = (0, 0, 1))
#
	#Test des partitions
    Partition_xm = object.cqshape.split(cq.Workplane().add(Plane_xm))
    Partition_xp = object.cqshape.split(cq.Workplane().add(Plane_xp))
    Partition_ym = object.cqshape.split(cq.Workplane().add(Plane_ym))
    Partition_yp = object.cqshape.split(cq.Workplane().add(Plane_yp))
    Partition_zm = object.cqshape.split(cq.Workplane().add(Plane_zm))
    Partition_zp = object.cqshape.split(cq.Workplane().add(Plane_zp))
    
    Partition=[Partition_xm.solids().all(), Partition_xp.solids().all(), Partition_ym.solids().all(), Partition_yp.solids().all(), Partition_zm.solids().all(), Partition_zp.solids().all()]

    for i in range (0,6):
        if len(Partition[i]) > 1:
            plist[i]=1
            object.add_pcount()

    for i in range(0,3):
        if plist[i*2] + plist[i*2+1] == 2:
            plist[i*2] = 0
            object.rem_pcount()
            plist[i*2+1] = 0
            object.rem_pcount()

    print(plist)

    if object.pcount==0: #inclusion ne dépasse pas du cube
        periodic_object.append(object.cqshape.solids().intersect(rve.Box))
    elif object.pcount==1:
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
    elif object.pcount==2:
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

    elif object.pcount==3: #inclusion sur un coin
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


    return periodic_object

def fuse_parts(solids, retain_edges):
            
    occ_Solids = solids[0].wrapped

    for i in range(1, len(solids)):
        fuse = BRepAlgoAPI_Fuse(occ_Solids, solids[i].wrapped)
        occ_Solids = fuse.Shape()

    if retain_edges == True:
        return cq.Solid(occ_Solids)
    else:
        upgrader = ShapeUpgrade_UnifySameDomain(occ_Solids, True, True, True)
        upgrader.Build();
        fixed = upgrader.Shape();

        return cq.Solid(fixed)

def cut_parts(solids):
    phase_cut=[]
    to_cuts=[]
    phase_cut.append(solids[-1])
    for i in range(len(solids)-1, 0, -1):
        occ_Solids = solids[i-1].wrapped
        to_cuts.append(solids[i])
        to_cut = fuse_parts(to_cuts, False)
        cut = BRepAlgoAPI_Cut(occ_Solids, to_cut.wrapped)
        phase_cut.append(cq.Solid(cut.Shape()))
    return phase_cut[::-1]
    
    

