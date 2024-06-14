# BUILD.PY
import os
import subprocess
import time
import sys
import numpy as np
import parmed as pmd

srcdir = "carlos_peptides/LONG/brush/"
Seq = {  # Molceule_Name : A_A_Sequence
    #'HHHG'  : 'ACE HIS HIS HIS GLY HIS HIS HIS GLY HIS HIS HIS GLY HIS HIS HIS GLY HIS HIS CHIS',
    #'HHHG+' : 'ACE HIP HIP HIP GLY HIP HIP HIP GLY HIP HIP HIP GLY HIP HIP HIP GLY HIP HIP CHIP',
    #'HHHHG' : 'ACE HIS HIS HIS HIS GLY HIS HIS HIS HIS GLY HIS HIS HIS HIS GLY HIS HIS HIS HIS GLY HIS HIS HIS CHIS',
    #'HHHHG+': 'ACE HIP HIP HIP HIP GLY HIP HIP HIP HIP GLY HIP HIP HIP HIP GLY HIP HIP HIP HIP GLY HIP HIP HIP CHIP',
    #'HHHH'  : 'ACE HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS CHIS',
    #'HHHH+' : 'ACE HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP CHIP',
    #'HHHP'  : 'ACE HIS HIS HIS PRO HIS HIS HIS PRO HIS HIS HIS PRO HIS HIS HIS PRO HIS HIS CHIS',
    #'HHHP+' : 'ACE HIP HIP HIP PRO HIP HIP HIP PRO HIP HIP HIP PRO HIP HIP HIP PRO HIP HIP CHIP',
    #'SSSS' : 'NSER SER SER SER SER SER SER SER SER SER NME',
    #'KEKE' : 'NLYS GLU LYS GLU LYS GLU LYS GLU LYS GLU NME',
    #'KEEK' : 'NLYS GLU GLU LYS GLU GLU LYS GLU GLU LYS NME'
    #'KKEE' : 'ACE GLY GLU GLU LYS LYS GLU GLU LYS LYS GLU GLU LYS LYS GLU GLU LYS LYS GLU GLU LYS LYS GLU GLU LYS LYS GLY NHE',
    #'KKKK' : 'ACE GLY GLU GLU GLU GLU LYS LYS LYS LYS GLU GLU GLU GLU LYS LYS LYS LYS GLU GLU GLU GLU LYS LYS LYS LYS GLY NHE',
    #    'EKGH' : 'NGLU LYS GLU LYS GLY HIS HIS HIS GLY HIS HIS HIS GLY HIS HIS CHIS',
    "EKEK": "NGLU LYS GLU LYS GLU LYS GLU LYS GLU CLYS",
}

for mol in Seq:
    print("----------  Starting Setup for {}  ----------\n".format(mol))
    # make local directory
    try:
        subprocess.call("mkdir -p {}{}".format(srcdir, mol), shell=True)
    except:
        print("Unable to make directory: {}".format(mol))
    try:
        subprocess.call(
            "cp Scripts/Others/sequence.tleap {}{}/".format(
                srcdir, mol
            ),  # LOCATION OF GLOBAL TLEAP FILES
            shell=True,
        )
    except:
        print("Unable to transfer tleap files")
    repLine = "struct = sequence {{{}}}".format(Seq[mol])
    # edit LEaP files
    with open("{}{}/sequence.tleap".format(srcdir, mol), "r") as file:
        data = file.readlines()
    # replace sequence LEaP file
    print("   Replacing line: {}\n    with: {}\n".format(data[5], repLine))
    data[5] = "{}\n".format(repLine)
    for i in range(len(data)):
        if "MOLECULE" in data[i]:
            data[i] = data[i].replace("MOLECULE", "{}{m}/{m}".format(srcdir, m=mol))
    with open("{}{}/sequence.tleap".format(srcdir, mol), "w") as file:
        file.writelines(data)
    print("LEaP Files Creation Successful\n")
    # Utilise TLEaP
    sub = subprocess.Popen(
        "tleap -f  {}{}/sequence.tleap".format(srcdir, mol),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
    )
    out, errors = sub.communicate()
    print(out.decode("ascii"))
    print("TLEaP Build run complete \n")
