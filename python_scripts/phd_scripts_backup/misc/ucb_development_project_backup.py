import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import os
import glob
import sys


srcdir = "./binding_trajectories"  # source directory
ref_traj = "5ai5"  # alignment | trajectory to align to
atm_ind = np.arange(240, 540)  # alignment | indices to align


# get list of trajectories to analyse
traj_list = [d.split("/")[-1] for d in glob.glob("{}/*".format(srcdir))]
print("Found trajectories: {}".format(traj_list))

# load the trajectories into dictionary, keys = directory names
trajectories = dict()
for t in traj_list:
    trajectories[t] = md.load_dcd(
        "{}/{}/traj.dcd".format(srcdir, t), top="{}/{}/input.pdb".format(srcdir, t)
    )
print("Successfully loaded : {}".format(trajectories.keys()))

# align the trajectories to the reference
for t in traj_list[1:]:
    try:
        trajectories[t].superpose(
            trajectories[traj_list[0]], frame=0, atom_indices=atm_ind, parallel=True
        )
    except:
        print("Alignment Error")
print("Trajectories successfully aligned.")

# define ligand centre of mass
