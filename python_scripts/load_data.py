d="""
===============================================================================
                                    LOADING DATA

===============================================================================
"""

import glob
import pickle
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def hills(filename):
    hills = [[], []]
    with open(filename) as f:
        lines = f.readlines()
        data = [l for l in lines if l[0] not in ("@", "#")]
        data = [[float(val) for val in line.split()] for line in data]
        hills[0] = ([l[0] for l in data])
        hills[1] = ([l[5] for l in data])
    return hills

def colvar(filename, output):
    #colvar_data = [[], [], []]
    with open(filename) as f:
        head = f.readlines()[0]
    head = head.split()[2:]
    # read in old COLVAR file into dataFrame
    # filters out comment lines and splits columns via whitespace
    old_col = pd.concat([df[df.time != "#!"] \
                        for df in pd.read_csv(filename,
                                              delim_whitespace=True,
                                              names=head,
                                              skiprows=1,
                                              chunksize=1000)])
    # round timestamps to ensure successful merging
    old_col['int_time'] = old_col['time'].astype(float).astype(int)
    # remove duplicate lines created by restarts
    old_col = old_col.drop_duplicates(subset='int_time', keep='last')
    if output == "as_pandas":
        return old_col
    elif output == "as_numpy":        
        return old_col.values.astype(float)

def fes(filename, is_rew):
    fes = [[], [], []]
    with open(filename) as f:
        lines = f.readlines()
        data = [l for l in lines if l[0] not in ("@", "#")]
    data = [[float(val) for val in line.split()] for line in data]
    breaks = [i for i, e in enumerate(data) if e == []]
    # remove blank lines
    for index in sorted(breaks, reverse=True):
        del data[index]
    # get number of bins and CV names from header
    if is_rew:
        nbins = int(breaks[0])
        [x_name,y_name] = ['RMSD to IN','RMSD to OUT']
    else: 
        nbins = int([l.split()[-1] for l in lines if "nbins_pp.proj" in l][0])
        [x_name, y_name] = lines[0].split()[2:4]
    # organise data into plotable arrays
    split_data = [data[i:i + nbins] for i in range(0, len(data), nbins)]
    z = []
    for block in split_data:
        z.append([l[2] for l in block])
        fes[1].append(block[0][1])
    fes[0] = [l[0] for l in split_data[0]]
    fes[2] = np.asarray(z)
    return fes, [x_name, y_name]

def fes_simple(filename, is_rew):
    fes = np.loadtxt(filename)
    return fes

def xvg(filename):
    """ load xvg """
    with open(filename) as f:
        lines = f.readlines()
        data = [l for l in lines if l[0] not in ("@", "#")]
    data = [[float(val) for val in line.split()] for line in data]
    return data

def cd(filename):
    """ load CD data """
    with open(filename) as f:
        lines = f.readlines()
        data = [l for l in lines if l[0] not in ("@", "#")]
    data = [[float(val) for val in line.split()] for line in data]
    return data
