from microgen.Operations import fuseParts
import cadquery as cq

def periodic(cqshape, rve):
    """ DESCRIPTION

    Parameters
    ----------
    cqshape : TYPE
        DESCRIPTION
    rve : TYPE
        DESCRIPTION

    Returns
    -------
    return_object_periodic[0].copy() : TYPE
        DESCRIPTION
    flat_list : TYPE
        DESCRIPTION
    """

    wk_plane = cq.Workplane().add(cqshape.Solids())
    periodic_object = []
    # plist = [0, 0, 0, 0, 0, 0]
    # pcount = 0
    face = ['x-', 'x+', 'y-', 'y+', 'z-', 'z+']
    intersected_faces = []

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

    # creation de faces de coupe
    Plane_xm = cq.Face.makePlane(basePnt=(0, 0, 0),      dir=(-1,  0,  0))
    Plane_xp = cq.Face.makePlane(basePnt=(rve.dx, 0, 0), dir=( 1,  0,  0))
    Plane_ym = cq.Face.makePlane(basePnt=(0, 0, 0),      dir=( 0, -1,  0))
    Plane_yp = cq.Face.makePlane(basePnt=(0, rve.dy, 0), dir=( 0,  1,  0))
    Plane_zm = cq.Face.makePlane(basePnt=(0, 0, 0),      dir=( 0,  0, -1))
    Plane_zp = cq.Face.makePlane(basePnt=(0, 0, rve.dz), dir=( 0,  0,  1))
    #
    # Test des partitions
    Partition_xm = wk_plane.split(cq.Workplane().add(Plane_xm))
    Partition_xp = wk_plane.split(cq.Workplane().add(Plane_xp))
    Partition_ym = wk_plane.split(cq.Workplane().add(Plane_ym))
    Partition_yp = wk_plane.split(cq.Workplane().add(Plane_yp))
    Partition_zm = wk_plane.split(cq.Workplane().add(Plane_zm))
    Partition_zp = wk_plane.split(cq.Workplane().add(Plane_zp))

    Partition = [
        Partition_xm.solids().all(),
        Partition_xp.solids().all(),
        Partition_ym.solids().all(),
        Partition_yp.solids().all(),
        Partition_zm.solids().all(),
        Partition_zp.solids().all(),
    ]

    for i in range(0, 6):
        if len(Partition[i]) > 1:
            intersected_faces.append(face[i])

    # A réfléchir
    if 'x-' in intersected_faces and 'x+' in intersected_faces:
        intersected_faces.remove('x-')
        intersected_faces.remove('x+')
    if 'y-' in intersected_faces and 'y+' in intersected_faces:
        intersected_faces.remove('y-')
        intersected_faces.remove('y+')
    if 'z-' in intersected_faces and 'z+' in intersected_faces:
        intersected_faces.remove('z-')
        intersected_faces.remove('z+')

    # print(intersected_faces)

    if len(intersected_faces) == 0:  # inclusion ne dépasse pas du cube
        periodic_object.append(wk_plane)
    elif len(intersected_faces) == 1: # intersection avec une face
        if intersected_faces == ['x-']:
            periodic_object.append(Partition_xm.solids(">X").intersect(rve.Box))
            periodic_object.append(
                Partition_xm.solids("<X").translate((rve.dx, 0, 0)).intersect(rve.Box)
            )
        elif intersected_faces == ['x+']:
            periodic_object.append(Partition_xp.solids("<X").intersect(rve.Box))
            periodic_object.append(
                Partition_xp.solids(">X").translate((-rve.dx, 0, 0)).intersect(rve.Box)
            )
        elif intersected_faces == ['y-']:
            periodic_object.append(Partition_ym.solids(">Y").intersect(rve.Box))
            periodic_object.append(
                Partition_ym.solids("<Y").translate((0, rve.dy, 0)).intersect(rve.Box)
            )
        elif intersected_faces == ['y+']:
            periodic_object.append(Partition_yp.solids("<Y").intersect(rve.Box))
            periodic_object.append(
                Partition_yp.solids(">Y").translate((0, -rve.dy, 0)).intersect(rve.Box)
            )
        elif intersected_faces == ['z-']:
            periodic_object.append(Partition_zm.solids(">Z").intersect(rve.Box))
            periodic_object.append(
                Partition_zm.solids("<Z").translate((0, 0, rve.dz)).intersect(rve.Box)
            )
        elif intersected_faces == ['z+']:
            periodic_object.append(Partition_zp.solids("<Z").intersect(rve.Box))
            periodic_object.append(
                Partition_zp.solids(">Z").translate((0, 0, -rve.dz)).intersect(rve.Box)
            )
    elif len(intersected_faces) == 2: # intersection with an edge

        if intersected_faces == ['x-', 'y-']:
            Partition_xpy = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
            )
            periodic_object.append(Partition_xpy.solids(">Y").intersect(rve.Box))
            periodic_object.append(
                Partition_xpy.solids("<Y").translate((0, rve.dy, 0)).intersect(rve.Box)
            )
            Partition_xmy = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
            )
            periodic_object.append(
                Partition_xmy.solids(">Y").translate((rve.dx, 0, 0)).intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmy.solids("<Y")
                .translate((rve.dx, rve.dy, 0))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x-', 'y+']:
            Partition_xpy = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_yp))
            )
            periodic_object.append(Partition_xpy.solids("<Y").intersect(rve.Box))
            periodic_object.append(
                Partition_xpy.solids(">Y").translate((0, -rve.dy, 0)).intersect(rve.Box)
            )
            Partition_xmy = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_yp))
            )
            periodic_object.append(
                Partition_xmy.solids("<Y").translate((rve.dx, 0, 0)).intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmy.solids(">Y")
                .translate((rve.dx, -rve.dy, 0))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x+', 'y-']:
            Partition_xmy = (
                cq.Workplane()
                .add(Partition_xp.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
            )
            periodic_object.append(Partition_xmy.solids(">Y").intersect(rve.Box))
            periodic_object.append(
                Partition_xmy.solids("<Y").translate((0, rve.dy, 0)).intersect(rve.Box)
            )
            Partition_xpy = (
                cq.Workplane()
                .add(Partition_xp.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
            )
            periodic_object.append(
                Partition_xpy.solids(">Y").translate((-rve.dx, 0, 0)).intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xpy.solids("<Y")
                .translate((-rve.dx, rve.dy, 0))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x+', 'y+']:
            Partition_xmy = (
                cq.Workplane()
                .add(Partition_xp.solids("<X"))
                .split(cq.Workplane().add(Plane_yp))
            )
            periodic_object.append(Partition_xmy.solids("<Y").intersect(rve.Box))
            periodic_object.append(
                Partition_xmy.solids(">Y").translate((0, -rve.dy, 0)).intersect(rve.Box)
            )
            Partition_xpy = (
                cq.Workplane()
                .add(Partition_xp.solids(">X"))
                .split(cq.Workplane().add(Plane_yp))
            )
            periodic_object.append(
                Partition_xpy.solids("<Y").translate((-rve.dx, 0, 0)).intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xpy.solids(">Y")
                .translate((-rve.dx, -rve.dy, 0))
                .intersect(rve.Box)
            )

        elif intersected_faces == ['x-', 'z-']:
            Partition_xpz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_zm))
            )
            periodic_object.append(Partition_xpz.solids(">Z").intersect(rve.Box))
            periodic_object.append(
                Partition_xpz.solids("<Z").translate((0, 0, rve.dz)).intersect(rve.Box)
            )
            Partition_xmz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_zm))
            )
            periodic_object.append(
                Partition_xmz.solids(">Z").translate((rve.dx, 0, 0)).intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmz.solids("<Z")
                .translate((rve.dx, 0, rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x-', 'z+']:
            Partition_xpz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_zp))
            )
            periodic_object.append(Partition_xpz.solids("<Z").intersect(rve.Box))
            periodic_object.append(
                Partition_xpz.solids(">Z").translate((0, 0, -rve.dz)).intersect(rve.Box)
            )
            Partition_xmz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_zp))
            )
            periodic_object.append(
                Partition_xmz.solids("<Z").translate((rve.dx, 0, 0)).intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmz.solids(">Z")
                .translate((rve.dx, 0, -rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x+', 'z-']:
            Partition_xmz = (
                cq.Workplane()
                .add(Partition_xp.solids("<X"))
                .split(cq.Workplane().add(Plane_zm))
            )
            periodic_object.append(Partition_xmz.solids(">Z").intersect(rve.Box))
            periodic_object.append(
                Partition_xmz.solids("<Z").translate((0, 0, rve.dz)).intersect(rve.Box)
            )
            Partition_xpz = (
                cq.Workplane()
                .add(Partition_xp.solids(">X"))
                .split(cq.Workplane().add(Plane_zm))
            )
            periodic_object.append(
                Partition_xpz.solids(">Z").translate((-rve.dx, 0, 0)).intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xpz.solids("<Z")
                .translate((-rve.dx, 0, rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x+', 'z+']:
            Partition_xmz = (
                cq.Workplane()
                .add(Partition_xp.solids("<X"))
                .split(cq.Workplane().add(Plane_zp))
            )
            periodic_object.append(Partition_xmz.solids("<Z").intersect(rve.Box))
            periodic_object.append(
                Partition_xmz.solids(">Z").translate((0, 0, -rve.dz)).intersect(rve.Box)
            )
            Partition_xpz = (
                cq.Workplane()
                .add(Partition_xp.solids(">X"))
                .split(cq.Workplane().add(Plane_zp))
            )
            periodic_object.append(
                Partition_xpz.solids("<Z").translate((-rve.dx, 0, 0)).intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xpz.solids(">Z")
                .translate((-rve.dx, 0, -rve.dz))
                .intersect(rve.Box)
            )

        elif intersected_faces == ['y-', 'z-']:
            Partition_ypz = (
                cq.Workplane()
                .add(Partition_ym.solids(">Y"))
                .split(cq.Workplane().add(Plane_zm))
            )
            periodic_object.append(Partition_ypz.solids(">Z").intersect(rve.Box))
            periodic_object.append(
                Partition_ypz.solids("<Z").translate((0, 0, rve.dz)).intersect(rve.Box)
            )
            Partition_ymz = (
                cq.Workplane()
                .add(Partition_ym.solids("<Y"))
                .split(cq.Workplane().add(Plane_zm))
            )
            periodic_object.append(
                Partition_ymz.solids(">Z").translate((0, rve.dy, 0)).intersect(rve.Box)
            )
            periodic_object.append(
                Partition_ymz.solids("<Z")
                .translate((0, rve.dy, rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['y-', 'z+']:
            Partition_ypz = (
                cq.Workplane()
                .add(Partition_ym.solids(">Y"))
                .split(cq.Workplane().add(Plane_zp))
            )
            periodic_object.append(Partition_ypz.solids("<Z").intersect(rve.Box))
            periodic_object.append(
                Partition_ypz.solids(">Z").translate((0, 0, -rve.dz)).intersect(rve.Box)
            )
            Partition_ymz = (
                cq.Workplane()
                .add(Partition_ym.solids("<Y"))
                .split(cq.Workplane().add(Plane_zp))
            )
            periodic_object.append(
                Partition_ymz.solids("<Z").translate((0, rve.dy, 0)).intersect(rve.Box)
            )
            periodic_object.append(
                Partition_ymz.solids(">Z")
                .translate((0, rve.dy, -rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['y+', 'z-']:
            Partition_ymz = (
                cq.Workplane()
                .add(Partition_yp.solids("<Y"))
                .split(cq.Workplane().add(Plane_zm))
            )
            periodic_object.append(Partition_ymz.solids(">Z").intersect(rve.Box))
            periodic_object.append(
                Partition_ymz.solids("<Z").translate((0, 0, rve.dz)).intersect(rve.Box)
            )
            Partition_ypz = (
                cq.Workplane()
                .add(Partition_yp.solids(">Y"))
                .split(cq.Workplane().add(Plane_zm))
            )
            periodic_object.append(
                Partition_ypz.solids(">Z").translate((0, -rve.dy, 0)).intersect(rve.Box)
            )
            periodic_object.append(
                Partition_ypz.solids("<Z")
                .translate((0, -rve.dy, rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['y+', 'z+']:
            Partition_ymz = (
                cq.Workplane()
                .add(Partition_yp.solids("<Y"))
                .split(cq.Workplane().add(Plane_zp))
            )
            periodic_object.append(Partition_ymz.solids("<Z").intersect(rve.Box))
            periodic_object.append(
                Partition_ymz.solids(">Z").translate((0, 0, -rve.dz)).intersect(rve.Box)
            )
            Partition_ypz = (
                cq.Workplane()
                .add(Partition_yp.solids(">Y"))
                .split(cq.Workplane().add(Plane_zp))
            )
            periodic_object.append(
                Partition_ypz.solids("<Z").translate((0, -rve.dy, 0)).intersect(rve.Box)
            )
            periodic_object.append(
                Partition_ypz.solids(">Z")
                .translate((0, -rve.dy, -rve.dz))
                .intersect(rve.Box)
            )

    elif len(intersected_faces) == 3:  # inclusion on a corner
        if intersected_faces == ['x-', 'y-', 'z-']:
            Partition_xpypz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(Partition_xpypz.solids(">Z").intersect(rve.Box))
            periodic_object.append(
                Partition_xpypz.solids("<Z")
                .translate((0, 0, rve.dz))
                .intersect(rve.Box)
            )
            Partition_xpymz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(
                Partition_xpymz.solids(">Z")
                .translate((0, rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xpymz.solids("<Z")
                .translate((0, rve.dy, rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmypz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(
                Partition_xmypz.solids(">Z")
                .translate((rve.dx, 0, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmypz.solids("<Z")
                .translate((rve.dx, 0, rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmymz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(
                Partition_xmymz.solids(">Z")
                .translate((rve.dx, rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmymz.solids("<Z")
                .translate((rve.dx, rve.dy, rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x-', 'y-', 'z+']:
            Partition_xpypz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(Partition_xpypz.solids("<Z").intersect(rve.Box))
            periodic_object.append(
                Partition_xpypz.solids(">Z")
                .translate((0, 0, -rve.dz))
                .intersect(rve.Box)
            )
            Partition_xpymz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(
                Partition_xpymz.solids("<Z")
                .translate((0, rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xpymz.solids(">Z")
                .translate((0, rve.dy, -rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmypz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(
                Partition_xmypz.solids("<Z")
                .translate((rve.dx, 0, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmypz.solids(">Z")
                .translate((rve.dx, 0, -rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmymz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(
                Partition_xmymz.solids("<Z")
                .translate((rve.dx, rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmymz.solids(">Z")
                .translate((rve.dx, rve.dy, -rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x-', 'y+', 'z-']:
            Partition_xpypz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(Partition_xpypz.solids(">Z").intersect(rve.Box))
            periodic_object.append(
                Partition_xpypz.solids("<Z")
                .translate((0, 0, rve.dz))
                .intersect(rve.Box)
            )
            Partition_xpymz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(
                Partition_xpymz.solids(">Z")
                .translate((0, -rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xpymz.solids("<Z")
                .translate((0, -rve.dy, rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmypz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(
                Partition_xmypz.solids(">Z")
                .translate((rve.dx, 0, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmypz.solids("<Z")
                .translate((rve.dx, 0, rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmymz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(
                Partition_xmymz.solids(">Z")
                .translate((rve.dx, -rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmymz.solids("<Z")
                .translate((rve.dx, -rve.dy, rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x-', 'y+', 'z+']:
            Partition_xpypz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(Partition_xpypz.solids("<Z").intersect(rve.Box))
            periodic_object.append(
                Partition_xpypz.solids(">Z")
                .translate((0, 0, -rve.dz))
                .intersect(rve.Box)
            )
            Partition_xpymz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(
                Partition_xpymz.solids("<Z")
                .translate((0, -rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xpymz.solids(">Z")
                .translate((0, -rve.dy, -rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmypz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(
                Partition_xmypz.solids("<Z")
                .translate((rve.dx, 0, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmypz.solids(">Z")
                .translate((rve.dx, 0, -rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmymz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(
                Partition_xmymz.solids("<Z")
                .translate((rve.dx, -rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmymz.solids(">Z")
                .translate((rve.dx, -rve.dy, -rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x+', 'y-', 'z-']:
            Partition_xpypz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(Partition_xpypz.solids(">Z").intersect(rve.Box))
            periodic_object.append(
                Partition_xpypz.solids("<Z")
                .translate((0, 0, rve.dz))
                .intersect(rve.Box)
            )
            Partition_xpymz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(
                Partition_xpymz.solids(">Z")
                .translate((0, rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xpymz.solids("<Z")
                .translate((0, rve.dy, rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmypz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(
                Partition_xmypz.solids(">Z")
                .translate((-rve.dx, 0, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmypz.solids("<Z")
                .translate((-rve.dx, 0, rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmymz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(
                Partition_xmymz.solids(">Z")
                .translate((-rve.dx, rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmymz.solids("<Z")
                .translate((-rve.dx, rve.dy, rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x+', 'y-', 'z+']:
            Partition_xpypz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(Partition_xpypz.solids("<Z").intersect(rve.Box))
            periodic_object.append(
                Partition_xpypz.solids(">Z")
                .translate((0, 0, -rve.dz))
                .intersect(rve.Box)
            )
            Partition_xpymz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(
                Partition_xpymz.solids("<Z")
                .translate((0, rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xpymz.solids(">Z")
                .translate((0, rve.dy, -rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmypz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(
                Partition_xmypz.solids("<Z")
                .translate((-rve.dx, 0, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmypz.solids(">Z")
                .translate((-rve.dx, 0, -rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmymz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(
                Partition_xmymz.solids("<Z")
                .translate((-rve.dx, rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmymz.solids(">Z")
                .translate((-rve.dx, rve.dy, -rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x+', 'y+', 'z-']:
            Partition_xpypz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(Partition_xpypz.solids(">Z").intersect(rve.Box))
            periodic_object.append(
                Partition_xpypz.solids("<Z")
                .translate((0, 0, rve.dz))
                .intersect(rve.Box)
            )
            Partition_xpymz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(
                Partition_xpymz.solids(">Z")
                .translate((0, -rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xpymz.solids("<Z")
                .translate((0, -rve.dy, rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmypz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(
                Partition_xmypz.solids(">Z")
                .translate((-rve.dx, 0, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmypz.solids("<Z")
                .translate((-rve.dx, 0, rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmymz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(
                Partition_xmymz.solids(">Z")
                .translate((-rve.dx, -rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmymz.solids("<Z")
                .translate((-rve.dx, -rve.dy, rve.dz))
                .intersect(rve.Box)
            )
        elif intersected_faces == ['x+', 'y+', 'z+']:
            Partition_xpypz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(Partition_xpypz.solids("<Z").intersect(rve.Box))
            periodic_object.append(
                Partition_xpypz.solids(">Z")
                .translate((0, 0, -rve.dz))
                .intersect(rve.Box)
            )
            Partition_xpymz = (
                cq.Workplane()
                .add(Partition_xm.solids("<X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(
                Partition_xpymz.solids("<Z")
                .translate((0, -rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xpymz.solids(">Z")
                .translate((0, -rve.dy, -rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmypz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids("<Y")
            )
            periodic_object.append(
                Partition_xmypz.solids("<Z")
                .translate((-rve.dx, 0, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmypz.solids(">Z")
                .translate((-rve.dx, 0, -rve.dz))
                .intersect(rve.Box)
            )
            Partition_xmymz = (
                cq.Workplane()
                .add(Partition_xm.solids(">X"))
                .split(cq.Workplane().add(Plane_ym))
                .solids(">Y")
            )
            periodic_object.append(
                Partition_xmymz.solids("<Z")
                .translate((-rve.dx, -rve.dy, 0))
                .intersect(rve.Box)
            )
            periodic_object.append(
                Partition_xmymz.solids(">Z")
                .translate((-rve.dx, -rve.dy, -rve.dz))
                .intersect(rve.Box)
            )

    occ_solids_list = [s.val().Solids() for s in periodic_object]
    flat_list = [item.copy() for sublist in occ_solids_list for item in sublist]
    to_fuse = [cq.Shape(s.wrapped) for s in flat_list]
    return_object_periodic = fuseParts(to_fuse, False)
    return (return_object_periodic[0].copy(), flat_list)