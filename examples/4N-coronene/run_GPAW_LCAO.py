from ase import *
from ase.visualize import *
from ase.io import*
from gpaw import *

import numpy as npy

mol = read('crazy_mol.traj')
# Molecule allready centered in the 15x15x10 cell
#cell = npy.loadtxt('input.lvs')
#mol.set_cell(cell)
#mol.set_pbc(False)
#mol.center()
xc='LDA'

view(mol)

calc = GPAW(txt='out_LCAO.txt',xc=xc,mode='lcao',basis='dzp')

mol.set_calculator(calc)

en = mol.get_potential_energy()

print en

calc.write('out_LCAO_'+xc+'.gpw',mode='all')

print
print "good bye"
print