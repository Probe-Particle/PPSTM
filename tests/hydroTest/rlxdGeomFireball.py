# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Pro-tip:
# This script will print out the relaxed structure from FIREBALL's "answer.xyz" file (that contains all of the relaxation steps) as a file called "relaxed.xyz". Enjoy!
# Run it with python3.x
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def tail(fname, n=1, bs=1024):
    f = open(fname)
    f.seek(0,2)
    l = 1-f.read(1).count('\n')
    B = f.tell()
    while n >= l and B > 0:
            block = min(bs, B)
            B -= block
            f.seek(B, 0)
            l += f.read(block).count('\n')
    f.seek(B, 0)
    l = min(l,n)
    lines = f.readlines()[-l:]
    f.close()
    return lines

# INPUT:
fname = "answer.xyz"

f = open(fname)
noAtoms = f.readline()
n = int(noAtoms)
print (n)
rlxd = tail(fname, n+2)

# OUTPUT:
with open('relaxed.xyz', 'w') as f:
	for line in rlxd:
		f.write(line)

