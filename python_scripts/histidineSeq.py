d="""
===============================================================================
                        HYDROPHOBIC SEQUENCE BUILDER
===============================================================================
"""

import numpy as np
import os
import subprocess
import glob
import sys
import argparse
import logging

########################################################
#                   PARSING INPUTS
########################################################

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=d, epilog=" ")

# required arguments
parser.add_argument("-x", type=str, default='GLY', help='Additional Inserted residue (default: %(default)s)')
# optional arguments
parser.add_argument("-r", type=int, default=3, help='Number of repitions in Hydrophilic block (default: %(default)s)', required=False)
parser.add_argument("-hisl", type=int, default=3, help='Histidine Lenght in block (default: %(default)s)',required=False)
# additional flags
parser.add_argument("-ncap", action="store_true", help='N Terminal cap applied? (default: %(default)s)')
parser.add_argument("-ccap", action="store_true", help='C Terminal cap applied? (default: %(default)s)')
parser.add_argument("-prot", action="store_true", help='PROTonation of histidines (default: %(default)s)')

# define variable from user input arguments
args = parser.parse_args()
xres    = args.x
reps    = args.r
his_len = args.hisl
his_prot = args.prot

########################################################
#               GENERATE A-A SEQUENCES
########################################################

seq =  (("HIS " if not his_prot else "HIP ")*his_len+ xres + " ") * reps
print(seq)
if args.ncap:
    seq = "ACE " + seq
else:
    seq = "N" + seq
print(seq)
if args.ccap:
    seq = seq + "NME"
else:
    seq =  seq[:-4] + "C" + seq[-4:]
print(seq)

#########################################################
#                   GENERATE TLEAP FILE
######################################################### 

fname = "seq_"+str(his_len)+"HIS+"+xres+"x"+str(reps)
if args.ncap:   fname = fname+"_Ncapped"
if args.ccap:   fname = fname+"_Ccapped"
if his_prot:    fname = fname+"_Prot"

ff_line = "source oldff/leaprc.ff14SB"
seq_line = "struct = sequence { " + seq + " }"
save_line= "savepdb struct " + fname + ".pdb\nquit"

with open('temp.tleap', 'w') as f:
    f.write(ff_line+"\n"+seq_line+"\n"+save_line)

print("Generated Tleap file. Will produce pbd: "+fname)

#########################################################
#                       RUN TLEAP
######################################################### 

try: 
    sub = subprocess.Popen( 'tleap -f ./temp.tleap',
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True 
    )
    out,errors =  sub.communicate() 
except:
    print("ERROR: TLeap unable to complete")
    



