import os
import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np
from time import sleep

#ib = 'rhys@ib-server.chem.ucl.ac.uk'
ib = 'rhys@login.archer.ac.uk'
ibdir = 'Main_Project/'
srcdir = 'carlos_peptides/short/'
mols = ['EKGH']
watermodel = [
        #"amber14sb_tip4pd.ff/"
        "amber_disp.ff/"
        #"charmm36m_sw.ff/"
        ]
ligand = False

#-----------------------------------------------------------------------------
folders = ['00-Prep','01-Min','02-NVT','03-NPT1','04-NPT2','05-MD']
for mol in mols:
    # make directory on IB Server
    try:
        subprocess.call(
            'ssh {} mkdir -p {}{}'.format(ib,ibdir,mol),
            shell=True)
    except:
        print( 'Unable to make directory: {}'.format(mol))
    try:
        subprocess.call(
            'scp -r Gromacs/* {}:~/{}{}/'.format(ib,ibdir,mol),
            shell=True)
    except:
        print( 'Unable to transfer files')    
    try:
        subprocess.call(
            'scp {}{m}/{m}_initial.pdb {}:~/{}{m}/00-Prep/{m}_initial.pdb'.format(srcdir,ib,ibdir,m=mol),
            shell=True) 
        if ligand: subprocess.call(
            'scp {}{m}/ligand*.itp {}:~/{}{m}/00-Prep/'.format(srcdir,ib,ibdir,m=mol),
            shell=True) 
        if ligand: subprocess.call(
            'scp {}{m}/FRG_GMX.gro {}:~/{}{m}/00-Prep/'.format(srcdir,ib,ibdir,m=mol),
            shell=True) 
        for f in folders: 
            subprocess.call(
                    'scp -r water_models/{} {}:~/{}{m}/{}'.format(watermodel[0],ib,ibdir,f,m=mol),
                shell=True) 
    except:
        print( 'Unable to transfer .gro and .top')
    print ('\nStarting batch submission... \n')
    # submit the new job
    sub = subprocess.Popen( 'ssh -T {}'.format(ib),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True, bufsize=0,
            shell=True )
    sub.stdin.write('cd {}{}/00-Prep/\n'.format(ibdir,mol))
    sub.stdin.write('pwd\n')
    #sub.stdin.write('qsub prep.sh\n')
    out, errors = sub.communicate()
    print (out)
    print (errors)
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

      


