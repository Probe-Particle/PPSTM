#!/env/python
from ase.io import read,write

ptcda=read('FFLJ_z.xsf')
write('input_plot.xyz',ptcda)
