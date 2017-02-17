from ase import *
from gpaw import *
from ase.io import *


import numpy as npy

xc='LDA'

slab, calc = restart('out_LCAO_'+xc+'.gpw')#,txt = 'dos_out.txt')
#e_fermi = calc.get_fermi_level()
#num_at = slab.get_number_of_atoms()
b=53
wf = calc.get_pseudo_wave_function(band=b)
fname=  'wf_'+str(b)+'-HOMO.cube'
print 'writing wf', b, 'to file', fname
write(fname, slab, data=wf)
b=54
wf = calc.get_pseudo_wave_function(band=b)
fname=  'wf_'+str(b)+'-LUMO.cube'
print 'writing wf', b, 'to file', fname
write(fname, slab, data=wf)


print
print "Bye bye"
print
