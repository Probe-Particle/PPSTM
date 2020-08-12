#!/usr/bin/python

import os
import numpy as np


# this library has functions for handling pre-STM analysis of an electronical state

# global variables:


# ==============================
# ============================== Pure python functions
# ==============================

def Lorentz(e,e0,eta):
    return eta/(2*np.pi* ((e-e0)**2 + 0.25*(eta**2)) );

def pPDOS(eig, coeffs, energies, eta=0.1, atoms=[], orbs='sp' ,spherical='all'):
    """
    pPDOS(eig, coeffs, energies, eta=0.1, atoms=[], orbs='sp' ,spherical='all')
    procedure for plotting pseudo-projected density of states; pseudo means that the density is not normalized, 
    since overlaps and of-diagonal terms are not treated in this procedure, and only valence electrons are read
    eig - eigen-numbers obtained from reading procedures
    coeffs - the LCAO coefficients obtained from the reading procedures
    energies - array of energies on which you want to calculate DOS - e.g energies = np.arange(-2.,2.,0.01)
    eta - width of the Lorentzian for smearing of the eigenstates
    atoms = [] = 'all'... all atoms; [0] 1st atom only; [1,5] 2nd & 6th atom ....
    orbs = 'sp' or 'spd'
    spherical = 'all' or 's' or 'p' or 'd' or 'px', 'py', 'pz', 'dxy', 'dxz', 'dyz', 'dz2', 'dx2y2'
    projection to only some of the spherical harmonics of the atomic orbitals
    """
    dim = coeffs.shape
    tmp= coeffs.flatten()
    sh = (int(dim[0]), int(dim[1]/4), 4) if (orbs=='sp') else (int(dim[0]), int(dim[1]/9), 9)
    coef = tmp.reshape(sh)
    if (orbs=='spd'):
        coef[:,:,4:9] *= 5 # rescale back what was lowered in the reading procedure
    PDOS = np.zeros(len(energies))

    if (spherical=='all'):
        sp = np.arange(4) if (orbs=='sp') else np.arange(9)
    elif (spherical=='s'):
        sp = np.array([0])
    elif (spherical=='p'):
        sp = np.arange(1,4)
    elif (spherical=='d'):
        sp = np.arange(4,9)
    elif (spherical=='py'):
        sp = np.array([1])
    elif (spherical=='pz'):
        sp = np.array([2])
    elif (spherical=='px'):
        sp = np.array([3])
    elif (spherical=='dxy'):
        sp = np.array([4])
    elif (spherical=='dyz'):
        sp = np.array([5])
    elif (spherical=='dz2'):
        sp = np.array([6])
    elif (spherical=='dxz'):
        sp = np.array([7])
    elif (spherical=='dx2y2'):
        sp = np.array([8])
    klist = np.arange(sh[1]) if (atoms==[] or atoms=='all') else np.array(atoms);
    #print "sp:",sp,"klist:",klist

    for i in range(len(energies)):
        lor = Lorentz(energies[i],eig[:],eta)
        for j in range(sh[0]):	# loop over eigen-states
            part_sum = 0;
            for k in klist:		# loop over atoms
                part_sum += np.sum(np.absolute(coef[j,k,sp])) #summing over different atomic orbitals
            PDOS[i] += lor[j]*part_sum
    
    #print PDOS
    return PDOS;


############## END OF LIBRARY ##################################
