#!/env/python
# Use this script to plot multiple ".dat" files in one figure. Enjoy!

from matplotlib import pyplot
from pylab import genfromtxt

mat0 = genfromtxt('new.kpts',skip_header=1)


fname = 'answer.xyz'

filein = open(fname)
nLines = len(list(filein))
print nLines
lines = filein.readlines()
nAtoms = lines[0].split()
nAtoms = [int(i) for i in nAtoms]
print nAtoms

with open(fname, 'r') as f:
	numLin = len(list(f))
	lines = f.read().splitlines()
	firstLine = lines[0]

    secondLine = lines[1]
    thirdLine = lines[2]
    fourthLine = lines[3]
    last_line = lines[-1]

    print firstLine
    print last_line

