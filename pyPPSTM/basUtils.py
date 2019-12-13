#!/usr/bin/python

import numpy as np
#import elements  
import math
from pylab import genfromtxt

# default variables:

default_atom_size     =  0.10

# procedures for loading geometry from different files:
name = 'answer.xyz'

def loadAtoms( name , sl=False):
	f = open(name,"r")
	n=0;
	l = f.readline()
	#print "--",l,"--"
	try:
		n=int(l)
	except:
                raise ValueError("First line of a xyz file should contain the "
                "number of atoms. Aborting...")
	if (n>0):
		n=int(l)
		e=[];x=[];y=[]; z=[]; q=[]
		i = 0;
		for line in f:
		    if sl :
			print " forced skipped line : ", line
			sl = False
		    else:
			words=line.split()
			nw = len( words)
			ie = None
			try:
				e.append( words[0] )
	   			x.append( float(words[1]) )
	   			y.append( float(words[2]) )
	   			z.append( float(words[3]) )
				if ( nw >=5 ):
					q.append( float(words[4]) )
				else:
					q.append( 0.0 )
			except:
				print " skipped line : ", line
	f.close()
	nDim = []
	lvec = [] 
	return [ e,x,y,z,q ], nDim, lvec

def loadGeometryIN( fname ):
	print "importin atoms from FHI-AIMS input"
	f = open(fname )
	e=[];x=[];y=[]; z=[]; q=[]
	lvec = [] 
	for i in range(10000):
		ws = f.readline().split()
		if (len(ws)>0):
			if (ws[0]=='atom'):
				e.append(ws[4]); x.append(float(ws[1])); y.append(float(ws[2])); z.append(float(ws[3])); q.append(0);
			elif (ws[0]=='lattice_vector'):
				lvec.append([float(ws[1]),float(ws[2]),float(ws[3])])
			elif (ws[0]=='trust_radius'):
				break
	f.close()
	#print "lvec", lvec
	#print "e,x,y,z", e,x,y,z
	nDim = []
	return [ e,x,y,z,q ], nDim, lvec

# other procedures for operating with geometries:

def multCell( xyz, cel, m=(2,2,1) ):
# For fireball theres a bug to be fixed:
# ETOT gets printed in the second line if you do relaxation, needs to be removed
# Otherwise it will be counted  as a part of es
	n = len(xyz[0])
	mtot = m[0]*m[1]*m[2]*n
	es = [None] * mtot
	xs = [None] * mtot
	ys = [None] * mtot
	zs = [None] * mtot
	j  = 0
	for ia in range(m[0]):
		for ib in range(m[1]):
			for ic in range(m[2]):
				dx = ia*cel[0][0] + ib*cel[1][0] + ic*cel[2][0]
				print "dx", dx
				dy = ia*cel[0][1] + ib*cel[1][1] + ic*cel[2][1]
				print "dy", dy
				dz = ia*cel[0][2] + ib*cel[1][2] + ic*cel[2][2]
				for i in xrange(n):
					es[j]=xyz[0][i]
					xs[j]=xyz[1][i] + dx
					ys[j]=xyz[2][i] + dy
					zs[j]=xyz[3][i] + dz
					j+=1
	return [es,xs,ys,zs]
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def loadKayz(name , sk):
	kayz = genfromtxt(name, skip_header= sk)
	print kayz, len(kayz), type(kayz), kayz.shape[1]
	if (len(kayz) > 0) and (type(kayz) == np.ndarray) and (kayz.shape[1] == 4):
		kx = kayz[:,0]
		ky = kayz[:,1]
		kz = kayz[:,2]
		w  = kayz[:,3]
		weights = np.array([kx,ky,kz])
		noWeights = np.array([kx,ky,kz])
		kayPoints = np.transpose(noWeights)
	return kayPoints, w




# =========== Utils for plotting atoms =========================

XSF_HEAD_0='''ATOMS
'''

XSF_HEAD_2='''
BEGIN_BLOCK_DATAGRID_3D                        
   some_datagrid      
   BEGIN_DATAGRID_3D_whatever 
'''

def At2XSF(atoms):
	XSF_HEAD_1=XSF_HEAD_0
	for i in range(len(atoms[0])):
		XSF_HEAD_1 = XSF_HEAD_1+str(atoms[0][i])+" "+str(atoms[1][i])+" "+str(atoms[2][i])+" "+str(atoms[3][i])+"\n "
	XSF_HEAD=XSF_HEAD_1 + XSF_HEAD_2
	#print "DEBUG: XSF_HEAD:"
	#print XSF_HEAD
	return XSF_HEAD ;
