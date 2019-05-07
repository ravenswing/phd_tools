d="""
===============================================================================
                        HYDROPHOBIC SEQUENCE BUILDER
===============================================================================
"""

import numpy as np
import os
import subprocess
import glob
import sys
import argparse
from itertools import product

########################################################
#                   PARSING INPUTS
########################################################

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=d, epilog=" ")

# required arguments
parser.add_argument("-pdb", type=str, default='mol.pdb', help='PDB file to arrange (default: %(default)s)')
parser.add_argument("-len", type=int, default=12, help='Number of residues in each peptide (default: %(default)s)')
parser.add_argument("-row", type=int, default=3, help='Number of rows of peptides (default: %(default)s)')
parser.add_argument("-col", type=int, default=3, help='Number of columns of peptides (default: %(default)s)')
parser.add_argument("-sep", type=float, default=0.0, help='Separation between each peptide "cylinder" - Angstroms (default: %(default)s)')
# optional arguments
parser.add_argument("-rad", type=float, default=7.0, help='Radius of the assumed peptide "cylinders" - Angstroms (default: %(default)s)')
#parser.add_argument("-hisl", type=int, default=3, help='Histidine Lenght in block (default: %(default)s)',required=False)
# additional flags
#parser.add_argument("-ncap", action="store_true", help='N Terminal cap applied? (default: %(default)s)')
#parser.add_argument("-ccap", action="store_true", help='C Terminal cap applied? (default: %(default)s)')
#parser.add_argument("-prot", action="store_true", help='PROTonation of histidines (default: %(default)s)')

# define variable from user input arguments
args = parser.parse_args()
pdb     = args.pdb
size    = [args.row,args.col]
sep     = args.sep

########################################################
#               COPY PBD ENOUGH TIMES
########################################################
try:
    fh = open(pdb, 'r')
except FileNotFoundError:
    print("ERROR: specified pdb not found")

stem = pdb[:-4]

try:
    subprocess.call("mkdir multi-pdb/",shell=True) 
except:
    print("ERROR: generating multi-pdbs directory")

for n in np.arange(np.prod(size)):
    try: 
        subprocess.call("cp {} multi-pdb/{}_{}.pdb".format(pdb,stem,n),shell=True)
    except:
        print("ERROR: unable to copy to multi-pdb directory")

#########################################################
#               GENERATE PYMOL CMD FILE
######################################################### 

line1 = "import pymol\nimport numpy as np\ncmd.delete('all')\n"
line2 = "for i in np.arange({}):\n    cmd.load('multi-pdb/{}_{{}}.pdb'.format(i))".format(np.prod(size), stem)
line3 = "\ncmd.orient('all')\n"

with open('temp.py','w') as f:
    f.write(line1)
    f.write(line2)
    f.write(line3)

#########################################################
#                       RUN TLEAP
######################################################### 
num = np.array(list(product([0,1,2], repeat=2)))
print(num)
num = np.stack(num)
num2 = np.sort(num,axis=1)
print(num2)

