#############################################
# 2019-05-29
# Crude script for formatting a PyMol selection to GROMACS index
# And printing it to stdout
# Giulio Mattedi
#############################################

import numpy as np
from pymol import cmd,stored

stored.list=[]
cmd.iterate("(sele)","stored.list.append(ID)")
idx = np.asarray(stored.list)

NUMS_PER_LINE = 15

print("[ index ]")
while True:
    print(
            ''.join(
                '%8d' % n for n in idx[:NUMS_PER_LINE]
                )
            )
    idx = idx[NUMS_PER_LINE:]

    if idx.size == 0: break