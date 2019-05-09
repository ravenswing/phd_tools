import subprocess
import argparse
import pandas as pd

import glob
import numpy as np

gmx_path = "/usr/local/gromacs/bin/gmx"

########################################################
#                   PARSING INPUTS
########################################################

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=d, epilog=" ")

# required arguments
parser.add_argument("-pdb", type=str, default='5ai0', help='PDB code (default: %(default)s)')
parser.add_argument("-fs", type=int, default=1, help='fs (default: %(default)s)')
parser.add_argument("-refpath", type=str, default='./reference.pdb', help='path to reference pdb for alignment in PyMol (default: %(default)s)')

# optional arguments
#parser.add_argument("-rad", type=float, default=7.0, help='Radius of the assumed peptide "cylinders" - Angstroms (default: %(default)s)')

# additional flags
#parser.add_argument("-prot", action="store_true", help='PROTonation of histidines (default: %(default)s)')

# define variable from user input arguments
args = parser.parse_args()
pdb     = args.pdb
fs      = args.fs

########################################################
#                  GENERATE OUT PDB
########################################################
wd = pdb+"_"+fs
stem = pdb+"-"fs

colvar = "{}/{}.COLVAR.old".format(wd,stem)
xtc = "{}/{}_rw_fin.xtc".format(wd,stem)
tpr = "{}/{}_ext.tpr".format(wd,pdb)

try:
    subprocess.call("cp ../UCB-350ns_Cterm/{}/{}.colvar {}".format(wd,pdb,colvar))
except:
    print("ERROR: cannot find old COLVAR.")

with open(colvar) as f: 
    lines = f.readlines()
    data = [l for l in lines if l[0] not in ("@", "#")]
    data = [l.split()[:2] for l in data]

df = pd.DataFrame({'time':[float(x[0]) for x in data[:]],'proj':[float(x[1]) for x in data[:]]})
df.sort_values('time',ascending=False, inplace=True)
df = df[ df['proj'].between(3.4,3.5,inclusive=True) ]

timestamp = round(df['time'].iloc[0],-1)
#print(timestamp)
#print(df.head(10))

out_name = "{}/{}_OUT.pdb".format(wd,stem)

try:
    subprocess.call("echo Protein_LIG | {} trjconv -s {} -f {} -o {} -b {t} -e {t} -n {}/i.ndx"\
            .format(gmx_path,tpr,xtc,out_name,wd,t=timestamp),shell=True)
except:
    print("ERROR: trjconv failed.")

########################################################
#               ALIGN ALL PDBS TO REF.
########################################################



