"""
===============================================================================
                                TRAJECTORY TOOLS
===============================================================================
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import parmed as pmd
import pytraj as pt
from glob import glob
from itertools import chain
import subprocess
import sys
import os


def _load_structure(in_str):
    # .pdb only req. 1 input
    if isinstance(in_str, str) and '.pdb' in in_str:
        return pt.load(in_str)
    # .r and same named .top
    elif isinstance(in_str, str):
        return pt.load(in_str, in_str.split('.')[0]+'.top')
    # explicitly specified crd and top
    elif isinstance(in_str, list) and len(in_str) == 2:
        return pt.load(in_str[0], in_str[1])
    # for references, can just be int i.e. frame no.
    elif isinstance(in_str, int):
        return in_str
    # if not any of the above raise an error
    else:
        raise ValueError("Structure not recognised")


def align(in_str, ref_str, out_str, aln_mask='@CA,C,N,O', strip_mask=None):
    # load the initial structure
    to_align = _load_structure(in_str)
    ref = _load_structure(ref_str)
    # run the alignment
    aligned = pt.align(to_align, mask=aln_mask, ref=ref)
    # if strip is required, perform the strip
    if strip_mask is not None:
        aligned = aligned.strip(strip_mask)
    # write the new str
    pt.write_traj(out_str, aligned, overwrite=True)


def cut_traj(trj_path, top, out_path, denom=100, split=False):
    full_trj = pt.iterload(trj_path, top)
    full_trj = full_trj.autoimage()
    print(f'Loaded trajectory: {trj_path}')
    if not split:
        start_point = 1
        print(f'NOT cutting traj so starting from {start_point}')
    else:
        start_point = int(full_trj.n_frames/2)
        print(f'CUTTING traj, starting from {start_point}')
    N = int(full_trj.n_frames/denom)
    print(f'Writing {N} frames')
    frames = np.linspace(start_point, full_trj.n_frames, num=N, dtype=int)-1
    pt.write_traj(out_path, full_trj, frame_indices=frames, overwrite=True)
    print(f'Saved new trajectory: {out_path}')


def measure_rmsd(trj_path, top_path, ref_str, rmsd_mask, aln_mask='@CA,C,N,O'):
    # load the trajectory w. topology
    traj = pt.iterload(trj_path, top_path)
    # load ref. structure if path is given
    ref = _load_structure(ref_str)
    # run autoimage to cluster and centre traj.
    traj = traj.autoimage()
    # align the traj. using backbone atoms
    traj = pt.align(traj, mask=aln_mask, ref=ref)
    # calculate rmsd
    data = pt.analysis.rmsd.rmsd_nofit(traj, mask=rmsd_mask, ref=ref)
    return data
