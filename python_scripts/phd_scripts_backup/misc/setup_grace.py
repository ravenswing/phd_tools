import os
import subprocess
import time
import sys


gr = "zcqsrev@grace.rc.ucl.ac.uk"
srcdir = "Hydrophobic/"
mols = ["HHHP", "HHHP+"]


for mol in mols:
    # make directory
    try:
        subprocess.call("ssh {} mkdir -p Scratch/{}".format(gr, mol), shell=True)
    except:
        print("Unable to make directory: {}".format(mol))
    # copy .gro and .top files to new directory
    try:
        subprocess.call(
            "scp {}{m}/05-MD/{m}_1barframe.gro {}:~/Scratch/{m}/".format(
                srcdir, gr, m=mol
            ),
            shell=True,
        )
        subprocess.call(
            "scp {}{m}/05-MD/{m}.top {}:~/Scratch/{m}/".format(srcdir, gr, m=mol),
            shell=True,
        )
    except:
        print("Unable to transfer top/gro files")
    # edit working directory in .sh run file
    with open("Scripts/Grace/mdrun.sh", "r") as file:
        data = file.readlines()
    #  print data
    print("Editing line:  {}".format(data[6]))
    data[1] = "#$ -N MD_{} \n".format(mol)
    data[6] = "#$ -wd /home/zcqsrev/Scratch/{} \n".format(mol)
    with open("Scripts/Grace/mdrun.sh", "w") as file:
        file.writelines(data)
    # edit working directory in .sh restart script
    with open("Scripts/Grace/mdrst.sh", "r") as file:
        data = file.readlines()
    print("Editing line:  {}".format(data[6]))
    data[6] = "#$ -wd /home/zcqsrev/Scratch/{} \n".format(mol)
    with open("Scripts/Grace/mdrst.sh", "w") as file:
        file.writelines(data)
    # copy (EDITED) .sh and .mdp files to Grace
    try:
        subprocess.call(
            "scp Scripts/Grace/md* {}:~/Scratch/{}/".format(gr, mol), shell=True
        )
    except:
        print("Unable to transfer mdrun files")
print("{} molecules successfully prepared".format(len(mols)))
# submit the first mdrun set to the Grace queue
print("\nStarting batch submission... \n")
for mol in mols:
    sub = subprocess.Popen(
        "ssh -T {}".format(gr),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        bufsize=0,
        shell=True,
    )

    sub.stdin.write("ls .\n")
    sub.stdin.write("cd Scratch/{}/\n".format(mol))
    sub.stdin.write("pwd\n")
    sub.stdin.write("qsub mdrun.sh\n")
    out, errors = sub.communicate()
    print(out)
    print(errors)
