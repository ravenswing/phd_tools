import os
import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np
from time import sleep


if int(sys.argv[1]) not in list(range(1,6)):
    print("\nERROR: Incorrect user input, please choose Stage 1-5\n")
    sys.exit()
else:
    stage = int(sys.argv[1])
    
#ibdir = 'baker_peptides/Amber-disp/' # MUST HAVE /
#ibdir = 'baker_peptides/CHARMM36m/'
ibdir = 'ucb_run/' 
mols = ['5ai5','5ai6','5akh','5alt','5am0']


ib = 'rhys@ib-server.chem.ucl.ac.uk'
proc = ['00-Prep', '01-Min','02-NVT','03-NPT1','04-NPT2','05-MD']
comm = {'00-Prep':'prep', '01-Min':'min','02-NVT':'nvt',
            '03-NPT1':'npt','04-NPT2':'npt2','05-MD':'md'}

for mol in mols:
    # make directory on IB Server 
    print ('\nStarting batch submission... \n')
    # submit the new job
    sub = subprocess.Popen( 'ssh -T {}'.format(ib),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True, bufsize=0,
            shell=True )
    sub.stdin.write('cd {}{}/{}/\n'.format(ibdir,mol,proc[stage]))
    sub.stdin.write('pwd\n')
    sub.stdin.write('module load gcc/7.2.0\nmodule load openmpi/gnu-7.2/2.1.2\n')
    sub.stdin.write('module load fftw/gnu-7.2/3.3.7\nmodule load gromacs\nmodule load plumed/gnu-7.2/2.4.1\n')
    sub.stdin.write('qsub {}.sh\n'.format(comm[proc[stage]]))
    out, errors = sub.communicate()
    print (out)
    print (errors) 

      


