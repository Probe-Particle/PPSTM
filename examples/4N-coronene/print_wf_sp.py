from ase import *
from gpaw import GPAW
import numpy as npy

# creates a Fireball-like PPSTM input files from GPAW runs

xc='LDA'
calc = GPAW('out_LCAO_'+xc+'.gpw',txt='tmp.txt')

slab = calc.get_atoms()
at_nos = slab.get_atomic_numbers()

n_at = slab.get_number_of_atoms()
n_bands = calc.get_number_of_bands()
pre_eig = calc.get_eigenvalues(kpt=0, spin=0, broadcast=True)

n_bands = 2

eig = pre_eig #- calc.get_fermi_level()
eig = eig[53:55]

with open('fermi.dat','w') as f:
    print ('Fermi Level =', calc.get_fermi_level() , 'eV', file=f)

with open('eigen.dat','w') as f:
    print ( "1" , n_bands, file=f)
    print ("  ---- The eigen values (Fermi = ", calc.get_fermi_level() ,"eV) ----", file=f)
    npy.savetxt(f,eig,fmt="%0.4f")

print ("n_at = ", n_at)
print ("n_bands = " , n_bands)
print ("eig:", eig)

# s-orb
X=npy.zeros((n_bands,n_at*2+1))
X[:,0] = eig

for i in range(53,55):
    h=0	
    for j in range(n_at):
        X[i-53,2*j+1] = calc.wfs.kpt_u[0].C_nM[i,h]
        h += calc.wfs.setups[j].nao

with open('phik_gpaw_s.dat','w') as f:
    print ( n_at, n_bands, calc.get_fermi_level(), file=f)
    npy.savetxt(f,X)

# py-orb
X=npy.zeros((n_bands,n_at*2+1))
X[:,0] = eig

for i in range(53,55):
    h=0	
    for j in range(n_at):
        if (at_nos[j] !=1) :		
            X[i-53,2*j+1] = calc.wfs.kpt_u[0].C_nM[i,h+1]
            h += calc.wfs.setups[j].nao

with open('phik_gpaw_py.dat','w') as f:
    print ( n_at, n_bands, calc.get_fermi_level(), file=f )
    npy.savetxt(f,X)

# pz-orb
X=npy.zeros((n_bands,n_at*2+1))
X[:,0] = eig

for i in range(53,55):
    h=0	
    for j in range(n_at):
        if (at_nos[j] !=1) :
            X[i-53,2*j+1] = calc.wfs.kpt_u[0].C_nM[i,h+2]
            h += calc.wfs.setups[j].nao

with open('phik_gpaw_pz.dat','w') as f:
    print ( n_at, n_bands, calc.get_fermi_level(), file=f )
    npy.savetxt(f,X)

# px-orb
X=npy.zeros((n_bands,n_at*2+1))
X[:,0] = eig

for i in range(53,55):
    h=0	
    for j in range(n_at):
        if (at_nos[j] !=1) :
            X[i-53,2*j+1] = calc.wfs.kpt_u[0].C_nM[i,h+3]
            h += calc.wfs.setups[j].nao

with open('phik_gpaw_px.dat','w') as f:
    print ( n_at, n_bands, calc.get_fermi_level(), file=f)
    npy.savetxt(f,X)

'''
#xy yz z2 xz x2-y2

# xy-orb
X=npy.zeros((n_bands,n_at*2+1))
X[:,0] = eig

f=open('phik_gpaw_dxy.dat','w')
print ( n_at, n_bands, file =f)
npy.savetxt(f,X,fmt="%4.2f")
f.close()

f=open('phik_gpaw_dyz.dat','w')
print ( n_at, n_bands, file =f)
npy.savetxt(f,X,fmt="%4.2f")
f.close()

f=open('phik_gpaw_dz2.dat','w')
print ( n_at, n_bands, file =f)
npy.savetxt(f,X,fmt="%4.2f")
f.close()


f=open('phik_gpaw_dxz.dat','w')
print ( n_at, n_bands, file =f)
npy.savetxt(f,X,fmt="%4.2f")
f.close()

f=open('phik_gpaw_dx2y2.dat','w')
print ( n_at, n_bands, file =f)
npy.savetxt(f,X,fmt="%4.2f")
f.close()
'''

print ("wave function written, keep on rocking")
