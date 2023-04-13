# Analyse the output of MD simulation using Driver
import os
import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np
from time import sleep

svr = 'pr1eeq09@mn3.bsc.es'
svrdir = '/gpfs/scratch/pr1eeq00/pr1eeq09/'
srcdir = 'Hydrophilic/'
mols = ['KEKE', 'KEEK','SSSS']
res = [1, 10]
# Available CVs = aRMSD
use = ['aRMSD']

#-----------------------------------------------------------------------------

CVs = {'aRMSD': 'a1: ALPHARMSD RESIDUES={}-{} TYPE=DRMSD LESS_THAN={{RATIONAL R_0=0.08 NN=8 MM=12}}'.format(res[0],res[1]),
        }

for mol in mols:

    try:
        subprocess.call(
        'scp Scripts/driver.sh {}:{}{m}/'.format(svr,svrdir,m=mol),
        shell=True)
        subprocess.call(
        'scp {}{m}/05-MD/{m}_protein.pdb {}:/gpfs/scratch/pr1eeq00/pr1eeq09/{m}/'.format(srcdir,svr,m=mol),
        shell=True)
    except:
        print( 'Unable to transfer driver.sh')
    try:
        print ('\nStarting batch submission... \n')
        # run trjconv on IB server          
        sub = subprocess.Popen( 'ssh -T {}'.format(svr),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True, bufsize=0,
            shell=True )
        sub.stdin.write('cd {}{}/\n'.format(svrdir,mol)) 
        sub.stdin.write('pwd\n')
        sub.stdin.write('sed -i -e \'s/MOLECULE/{}/\' driver.sh\n'.format(mol))
        sub.stdin.write('rm -f -- plumed.dat\n')
        sub.stdin.write('echo \'MOLINFO STRUCTURE={m}_protein.pdb\n\' >> plumed.dat\n'.format(m=mol))
        for cv in CVs:
            if cv not in use: continue
            sub.stdin.write('echo \'{}\n\' >> plumed.dat\n'.format(CVs[cv]))
              
        sub.stdin.write('echo \'PRINT ARG=* STRIDE=1 FILE=COLVAR\n\' >> plumed.dat\n')
        sub.stdin.write('bash driver.sh\n')
        out, errors = sub.communicate()
        print (out)
        print (errors)
    except:
        print( 'Error in driver setup/run.')
    print( 'Successfully run Driver for {}  |  Attempting rsync...'.format(mol))
    try:
        subprocess.call(
        'rsync -rvzhPe ssh --exclude "*.trr" --ignore-existing {}:{}{m}/ {}{m}/06-RunMD/'.format(svr,svrdir,srcdir,m=mol),
        shell=True)
    except:
        print( 'Error in RSYNC')
             

    
'''
wait = True
finished = []
while wait:
    time.sleep(60)
    for mol in mols:
        if mol in finished: continue
        ls = subprocess.Popen(['ssh','{}'.format(ib), 'ls {}/01-Min/'.format(mol)],
             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err =  ls.communicate()
        output = out.decode('utf-8').splitlines()
        if '{}_min_energy.xvg'.format(mol) in output:
            try:
                subprocess.call(
                    'scp {}{m}/{m}.gro {}:~/{m}/01-Min/'.format(srcdir,ib,m=mol),
                    shell=True)
            except:
                print( 'Unable to transfer BACK .xvg file') 
        if len(finished) == len(mols): wait = False    

    '''

      


