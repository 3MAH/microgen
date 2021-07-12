import os
import numpy as np

from microgen.Functions import *
from microgen.Sphere import *
from microgen.Cylinder import *
from microgen.Ellipsoid import *
from microgen.Bar import *
from microgen.Rve import *

class phase:

    def __init__(self,number,shape,xc,yc,zc,psi,theta,phi,param_geom,path_data):
        self.number=number
        self.shape=shape
        self.xc=xc
        self.yc=yc
        self.zc=zc
        self.psi=psi
        self.theta=theta
        self.phi=phi
        self.param_geom=param_geom
        self.path_data=path_data
        
        self.center=np.array([self.xc,self.yc,self.zc])
        self.angle=np.array([self.psi,self.theta,self.phi])
        self.name=self.shape+str(self.number)
        self.pcount=0

    def __cmp__(self, other):
        return cmp(self.number, other.number)

    def add_pcount(self):
        self.pcount = self.pcount+1

    def rem_pcount(self):
        self.pcount = self.pcount-1
    
#----------GENERATE PHASES----------------------------------------------------------------------------------

    def generate(self,rve):

        if self.shape == 'cylinder' :
            objet = cylinder(self.center,self.angle,self.param_geom[0],self.param_geom[1],self.number)
            objet = objet.create_cylinder()
        if self.shape == 'bar' :
            objet = bar(self.center,self.angle,self.param_geom[0],self.param_geom[1],self.number)
            objet = objet.create_bar()
        if self.shape == 'sphere' :
            objet = sphere(self.center,self.param_geom[0],self.number)
            objet = objet.create_sphere()
        if self.shape == 'ellipsoid' :
            objet = ellipsoid(self.center,self.angle,self.param_geom[0],self.param_geom[1],self.param_geom[2],self.number)
            objet = objet.create_ellipsoid()
#        if self.shape == 'tpms' :
#            objet = tpms(self.center,self.param_geom[0],self.param_geom[1],self.param_geom[2],self.number)
#            objet = objet.create_tpms(self.path_data,rve)

        self.cqshape = objet

def read_phases(path_data,phases_file,phases):

    nphases = 0
    cnt_phase = 0
    nprops = []
    buf = ""
    path_inputfile = path_data + phases_file
    remove_empty_lines(path_inputfile)
    
    try:
        fp = open(path_inputfile)
        for line in enumerate(fp):
            if line == '\n':
                nphases -= 1
            else:
                row_split = line[1].split()
                if nphases > 0:
                    
                    nprops_phase = 0
                    if(row_split[1] == 'matrix'):
                        nprops_phase = 0
                    elif(row_split[1] == 'sphere'):
                        nprops_phase = 1
                    elif(row_split[1] == 'cylinder'):
                        nprops_phase = 2
                    elif(row_split[1] == 'bar'):
                        nprops_phase = 2
                    elif(row_split[1] == 'ellipsoid'):
                        nprops_phase = 3
                    elif(row_split[1] == 'tpms'):
                        nprops_phase = 3
                                                                    
                    nprops.append(nprops_phase)
                nphases += 1
    finally:
        fp.seek(0)
        fp.close()

    print(path_data)
    print(phases_file)
    

    try:
        fp = open(path_inputfile)
        for line in enumerate(fp):
            print (line)
            if line == '\n':
                cnt_phase -= 1
            else:
                row_split = line[1].split()
                if cnt_phase > 0:
                    
                    number = int(row_split[0])
                    shape = row_split[1]
                    xc = float(row_split[2])
                    yc = float(row_split[3])
                    zc = float(row_split[4])
                    psi = float(row_split[5])
                    theta = float(row_split[6])
                    phi = float(row_split[7])

                    props = []
                    if shape == 'tpms':
                        for prop in range(0,nprops[cnt_phase-1]):
                            props.append(row_split[8+prop])
                    else:
                        for prop in range(0,nprops[cnt_phase-1]):
                            props.append(float(row_split[8+prop]))
                                                
                    pha = phase(number,shape,xc,yc,zc,psi,theta,phi,props,path_data)
                    phases.append(pha)
                    print (phases[-1].shape)
                cnt_phase += 1
    finally:
        fp.close()
        
