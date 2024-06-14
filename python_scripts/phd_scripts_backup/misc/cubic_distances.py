# Cubic

import mdtraj as md
import numpy as np
import pandas as pd
from itertools import product


def getMating(oxy_atoms, nit_atoms, num_pep):
    # Add atoms numbers to dictionary
    oxygens = {}
    nitrogens = {}
    for i in range(num_pep):
        oxygens[i] = oxy_atoms[i * 10 : i * 10 + 10]
        nitrogens[i] = nit_atoms[i * 5 : i * 5 + 5]

    atom_pairs = np.array([])
    this_pair = np.array([])
    expected = 0

    # For each peptide-group of nitrogens...
    for n in nitrogens:
        # Make a list of peptide-groups *excluding* the current group
        partners = list(nitrogens.keys())
        partners.remove(n)

        # For each peptide-group of oxygens in that list...
        for o in partners:
            # Load in the relevant atom numbers from the dictionary
            oxy_half = oxygens[o]
            nit_half = nitrogens[n]

            # Create a product of the arrays, containing every combination
            this_pair = np.array(list(product(oxy_half, nit_half)))

            # Add the product to the total list of products
            if len(atom_pairs) == 0:
                atom_pairs = this_pair
            else:
                atom_pairs = np.concatenate((atom_pairs, this_pair))

    return atom_pairs


srcdir = "./"

colnames = ["residue", "atom"]
# Load in oxygen data
oxy_list = pd.read_csv(
    "{}/atom_refs/cubic_oxy.txt".format(srcdir), names=colnames, delimiter="\t"
)
oxy_residues = np.array(oxy_list.residue.tolist())
oxy_atoms = np.array(oxy_list.atom.tolist())
ordered_oxy_atoms = np.linspace(1, len(oxy_atoms), len(oxy_atoms))

# Load in nitrogen data
nit_list = pd.read_csv(
    "{}/atom_refs/cubic_nit.txt".format(srcdir), names=colnames, delimiter="\t"
)
nit_residues = np.array(nit_list.residue.tolist())
nit_atoms = np.array(nit_list.atom.tolist())
ordered_nit_atoms = np.linspace(1, len(nit_atoms), len(nit_atoms))

atom_pairs = getMating(oxy_atoms, nit_atoms, 9)
atom_order = getMating(ordered_oxy_atoms, ordered_nit_atoms, 9)

# iad = ['05a-MD','05b-MD','xyz','z']
iad = ["05a-MD", "05b-MD"]
segments = np.linspace(0, 450, 10)

for r in range(len(iad)):
    mol = iad[r]

    for i in range(len(segments)):
        lower = int(segments[i]) + 1
        upper = int(segments[i]) + 50

        # Compute contacts for the given trajectory
        target = "cut_{l}-{u}ns.xtc".format(l=lower, u=upper)
        traj = md.load_xtc(
            "{}/{m}/{t}".format(srcdir, m=mol, t=target),
            top="{}/{m}/{m}_protein.gro".format(srcdir, m=mol),
        )

        # contacts = md.compute_contacts(traj, scheme='sidechain')
        interactions = md.compute_distances(traj, atom_pairs, periodic=False)

        # Save plottable variables
        np.save(
            "{}/{m}/{m}_{l}-{u}_interactions.npy".format(
                srcdir, m=mol, l=lower, u=upper
            ),
            interactions,
        )

    # np.save('{}/{m}/{m}_atom_pairs.npy'.format(srcdir, m=mol), atom_pairs)
    np.save("{}/{m}/{m}_atom_order.npy".format(srcdir, m=mol), atom_order)
