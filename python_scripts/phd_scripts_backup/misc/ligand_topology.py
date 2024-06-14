import subprocess

srcdir = "ucb_development_project/BindingTrajectory/"
mols = ["5ai6", "5alt", "5akh", "5am0"]
filename = "FRG_GMX.top"

for mol in mols:
    with open("{}{}/{}".format(srcdir, mol, filename), "r") as f:
        d = "["
        s = [d + e for e in f.read().split(d) if e]

    # put atomtypes in separate itp file
    with open("{}{}/ligand_atomtypes.itp".format(srcdir, mol), "w+") as f:
        f.write([st for st in s if "atomtypes" in st][0])

    # remove all necessary sections from ligand topology
    remove = ["defaults", "atomtypes", "molecules", "system"]
    include = "\n".join([st for st in s[1:] if not any(x in st for x in remove)])
    with open("{}{}/ligand.itp".format(srcdir, mol), "w+") as f:
        f.write(include)
    print("  {}  |  Ligand topology successfully split".format(mol))
