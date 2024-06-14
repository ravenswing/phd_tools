# BUILD.PY
import os
import subprocess
import time
import sys
import numpy as np
import parmed as pmd


ib = "rhys@ib-server.chem.ucl.ac.uk"
srcdir = "Hydrophobic_Capped/"
Seq = {  # Molceule_Name : A_A_Sequence
    #'HHHG'  : 'ACE HIS HIS HIS GLY HIS HIS HIS GLY HIS HIS HIS GLY HIS HIS HIS GLY HIS HIS CHIS',
    #'HHHG+' : 'ACE HIP HIP HIP GLY HIP HIP HIP GLY HIP HIP HIP GLY HIP HIP HIP GLY HIP HIP CHIP',
    #'HHHH'  : 'ACE HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS HIS CHIS',
    #'HHHH+' : 'ACE HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP HIP CHIP',
    #'HHHP'  : 'ACE HIS HIS HIS PRO HIS HIS HIS PRO HIS HIS HIS PRO HIS HIS HIS PRO HIS HIS CHIS',
    "HHHP+": "ACE HIP HIP HIP PRO HIP HIP HIP PRO HIP HIP HIP PRO HIP HIP HIP PRO HIP HIP CHIP",
    #'SSSS' : 'NSER SER SER SER SER SER SER SER SER CSER',
    #'KEKE' : 'NLYS GLU LYS GLU LYS GLU LYS GLU LYS CGLU',
    #'KEEK' : 'NLYS GLU GLU LYS GLU GLU LYS GLU GLU CLYS'
}
conc = 0.150  # molar : concentration


for mol in Seq:
    print("----------  Starting Setup for {}  ----------\n".format(mol))
    # make local directory
    try:
        subprocess.call("mkdir -p {}{}".format(srcdir, mol), shell=True)
    except:
        print("Unable to make directory: {}".format(mol))
    try:
        subprocess.call(
            "cp Scripts/build.tleap {}{}/".format(
                srcdir, mol
            ),  # LOCATION OF GLOBAL TLEAP FILES
            shell=True,
        )
        subprocess.call("cp Scripts/solv.tleap {}{}/".format(srcdir, mol), shell=True)
    except:
        print("Unable to transfer tleap files")
    repLine = "struct = sequence {{ {} }}".format(Seq[mol])
    # repLine = 'struct = loadpdb IB_output/{}_compact.pdb'.format(mol)
    # edit LEaP files
    for fl in ["build", "solv"]:
        with open("{}{}/{}.tleap".format(srcdir, mol, fl), "r") as file:
            data = file.readlines()
        # replace sequence LEaP file
        print(
            "   Replacing line in {}: {}\n    with: {}\n".format(fl, data[5], repLine)
        )
        data[5] = "{}\n".format(repLine)
        for i in range(len(data)):
            if "MOLECULE" in data[i]:
                data[i] = data[i].replace("MOLECULE", "{}{m}/{m}".format(srcdir, m=mol))
        with open("{}{}/{}.tleap".format(srcdir, mol, fl), "w") as file:
            file.writelines(data)
    print("LEaP Files Creation Successful\n")
    # Solvate to get volume of box...
    sub = subprocess.Popen(
        "tleap -f  {}/{}/solv.tleap".format(srcdir, mol),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
    )
    out, errors = sub.communicate()
    output = out.decode("utf-8").splitlines()
    for line in output:
        print(line)
        if "Volume" in line:
            V = float(line.split(" ")[3])
    print("Solvated with volume: {} \n".format(V))
    N = round(conc * (V * 1e-27) * 6.023e23)
    print(N)
    # add NaCl ions to reach desired concentration
    with open("{}{}/build.tleap".format(srcdir, mol), "r") as file:
        data = file.readlines()
    data[11] = "addIons2 struct Cl- {}\n".format(N)
    data[12] = "addIons2 struct Na+ {}\n".format(N)
    with open("{}{}/build.tleap".format(srcdir, mol), "w") as file:
        file.writelines(data)
    print("Added {} NaCl ions to reach {} moldm3 concentration.\n".format(N, conc))
    # Utilise TLEaP
    sub = subprocess.Popen(
        "tleap -f  {}/{}/build.tleap".format(srcdir, mol),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
    )
    out, errors = sub.communicate()
    print(out)
    print(errors)
    print("TLEaP Build run complete, starting conversion \n")
    # Amber -> Gromace Conversion
    traj = pmd.load_file(
        "{}{m}/{m}.prmtop".format(srcdir, m=mol), "{}{m}/{m}.rst7".format(srcdir, m=mol)
    )
    traj.save("{}{m}/{m}.gro".format(srcdir, m=mol))
    traj.save("{}{m}/{m}.top".format(srcdir, m=mol))
    print(
        "Amber -> Gromacs Conversion Complete\n \n Setup Completed for {}\n".format(mol)
    )
