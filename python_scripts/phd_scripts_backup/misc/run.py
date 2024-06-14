import os
import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np
from time import sleep

ib = "rhys@ib-server.chem.ucl.ac.uk"
ibdir = "TIP4P-Ew/"
srcdir = "Hydrophilic_TIP4PEw/"
mols = ["KKKK", "KKEE"]
tol = 0.1  # percentage tolerance for min check:


def isMin(x, y, tol, xvgfile):
    ymax = max(y)
    endi = len(data) - int(np.floor(len(data) / 10))
    final10 = y[endi:]
    std = np.std(final10)
    if std < np.absolute(ymax * tol * 0.01):
        return True
    else:
        return False
    with open(xvgfile) as f:
        lines = f.readlines()
    data = [l for l in lines if l[0] not in ("@", "#")]
    data = [[float(val) for val in line.split()[:2]] for line in data]
    x, y = [l[0] for l in data], [l[1] for l in data]
    print(len(data))


# -----------------------------------------------------------------------------

for mol in mols:
    # make directory on IB Server
    try:
        subprocess.call("ssh {} mkdir -p {}{}".format(ib, ibdir, mol), shell=True)
    except:
        print("Unable to make directory: {}".format(mol))
    try:
        subprocess.call(
            "scp -r Gromacs/* {}:~/{}{}/".format(ib, ibdir, mol), shell=True
        )
    except:
        print("Unable to transfer files")
    try:
        subprocess.call(
            "scp {}{m}/00-Prep/{m}_m4.gro {}:~/{}{m}/01-Min/{m}.gro".format(
                srcdir, ib, ibdir, m=mol
            ),
            shell=True,
        )
        subprocess.call(
            "scp {}{m}/00-Prep/{m}_m4.top {}:~/{}{m}/01-Min/{m}.top".format(
                srcdir, ib, ibdir, m=mol
            ),
            shell=True,
        )
    except:
        print("Unable to transfer .gro and .top")
    print("\nStarting batch submission... \n")
    # submit the new job
    sub = subprocess.Popen(
        "ssh -T {}".format(ib),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        bufsize=0,
        shell=True,
    )
    sub.stdin.write("cd {}{}/01-Min/\n".format(ibdir, mol))
    sub.stdin.write("pwd\n")
    sub.stdin.write("qsub min.sh\n")
    out, errors = sub.communicate()
    print(out)
    print(errors)
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
