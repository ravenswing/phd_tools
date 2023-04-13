import os
import subprocess
import time
import sys
import parmed as pmd

gr = 'zcqsrev@grace.rc.ucl.ac.uk'

mols = ['HHHH', 'HHHH+']

for mol in mols:
    #subprocess.call( 'ssh {} cd Scratch/{}/ && qsub mdrun.sh'.format(gr,mol), shell=True)
    sub = subprocess.Popen( 'ssh -T {}'.format(gr),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        universal_newlines=True, bufsize=0,
        shell=True )

    sub.stdin.write('ls .\n')

    sub.stdin.write('cd Scratch/{}/\n'.format(mol))

    sub.stdin.write('pwd\n')

    sub.stdin.write('qsub mdrun.sh\n')
    out, errors = sub.communicate()

    print (out)
    print (errors)
 
