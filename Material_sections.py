# -*- coding: utf8 -*-
# Do not delete the following import lines
import numpy as np

class MatSection:

	def __init__(self,number,name,umat_name,psi_mat,theta_mat,phi_mat,nprops,nstatev,props):
		self.number=number	#Number of the section
		self.name=name		# Name of the section
		self.umat_name=umat_name	# Name of the constitutive law
		self.psi_mat=psi_mat		# Name of the Euler angle for rotation around the z axis
		self.theta_mat=theta_mat	# Name of the Euler angle for rotation around the x' axis
		self.phi_mat=phi_mat		# Name of the Euler angle for rotation around the z'' axis
		self.nprops=nprops			# Number of material parameters
		self.nstatev=nstatev		# Number of internal variables for user constitutive laws
		self.props=props			# Table of the material properties
			
		c1 = np.cos(self.psi_mat*np.pi/180.0)
		s1 = np.sin(self.psi_mat*np.pi/180.0)
		c2 = np.cos(self.theta_mat*np.pi/180.0)
		s2 = np.sin(self.theta_mat*np.pi/180.0)
		c3 = np.cos(self.phi_mat*np.pi/180.0)				
		s3 = np.sin(self.phi_mat*np.pi/180.0)
		self.R = np.array([[c3*c2*c1-s3*s1,-c1*s3-c3*c2*s1,c3*s2], [c3*s1+c2*c1*s3,c3*c1-c2*s3*s1,s3*s2], [-c1*s2,s2*s1,c2]])
	
		x = np.array([[1.0],[0.0],[0.0]])	
		y = np.array([[0.0],[1.0],[0.0]])
		x2mat = np.ndarray.flatten(self.R.dot(x))
		y2mat = np.ndarray.flatten(self.R.dot(y))
		
		self.p0 = (0.0, 0.0, 0.0) 
		self.p1 = x2mat.tolist()
		self.p2 = y2mat.tolist()			
	
