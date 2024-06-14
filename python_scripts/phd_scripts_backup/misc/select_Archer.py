import os
import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np
from time import sleep

ib = "rhys@login.archer.ac.uk"
ibdir = "/work/e280/e280/rhys/Main_Project/"
srcdir = "carlos_peptides/short/"
mols = ["EKGH"]

# -----------------------------------------------------------------------------

for mol in mols:
    # load the data from pressure data file
    data = np.loadtxt("{}{m}/05-MD/{m}_pressure1bar.dat".format(srcdir, m=mol))
    ind = 0
    dif = 0.5
    # select the timestamp of the frame with pressure closest to 1 bar
    for i in range(len(data)):
        d = np.absolute(1.0 - data[i][1])
        if d < dif:
            ind = i
            dif = d
        else:
            continue
    ts = int(data[ind][0])
    print(
        "\n Selected Timestamp:  {}         |  Minimum Difference:  {}".format(ts, dif)
    )
    print("\nStarting batch submission... \n")
    # run trjconv on IB server
    sub = subprocess.Popen(
        "ssh -T {}".format(ib),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        bufsize=0,
        shell=True,
    )
    sub.stdin.write("cd {}{}/05-MD/\n".format(ibdir, mol))
    sub.stdin.write("pwd\n")
    sub.stdin.write(
        "echo '#!/bin/bash\n#$ -cwd -V\n#$ -q all.q\n#$ -pe openmpi 1\n' >> select.sh\n"
    )
    sub.stdin.write(
        "echo 'export LD_LIBRARY_PATH=/usr/local/plumed-2.3.0/lib:$LD_LIBRARY_PATH\n' >> select.sh\n"
    )
    sub.stdin.write(
        "echo 'export GMX=/usr/local/gromacs/5.1.4-plumed-2.3.0/bin/gmx_mpi\n' >> select.sh\n"
    )
    sub.stdin.write(
        "echo 'echo 0 | $GMX trjconv -s md.tpr -f md.trr -o {m}_1barframe.gro -b {t} -e {t}' >> select.sh\n".format(
            m=mol, t=ts
        )
    )
    sub.stdin.write("qsub select.sh\n")
    out, errors = sub.communicate()
    print(out)
    print(errors)
    #

"""
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

    """
