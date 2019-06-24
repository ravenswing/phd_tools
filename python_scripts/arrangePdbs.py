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
import pandas as pd
from itertools import product

########################################################
#                   PARSING INPUTS
########################################################

parser = argparse.ArgumentParser(\
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description=d, epilog=" ")

# required arguments
parser.add_argument("-pdb", type=str, default='mol.pdb',
                    help='PDB file to arrange (default: %(default)s)')
parser.add_argument("-len", type=int, default=12,
                    help='Number of residues in each peptide '+
                        '(default: %(default)s)')
parser.add_argument("-row", type=int, default=3,
                    help='Number of rows of peptides (default: %(default)s)')
parser.add_argument("-col", type=int, default=3,
                    help='Number of columns of peptides (default: %(default)s)')
parser.add_argument("-sep", type=float, default=0.0,
                    help='Separation between each peptide "cylinder" '+
                        '- Angstroms (default: %(default)s)')
# optional arguments
parser.add_argument("-rad", type=float, default=7.0,
                    help='Radius of the assumed peptide "cylinders" '+
                        '- Angstroms (default: %(default)s)')
# define variable from user input arguments
args = parser.parse_args()
pdb     = args.pdb
size    = [args.row,args.col]

########################################################
#               COPY PBD ENOUGH TIMES
########################################################

# confirm that the target pdb exists
try:
    fh = open(pdb, 'r')
except FileNotFoundError:
    print("ERROR: specified pdb not found")
# stem = filename - used for all subsequent files
stem = pdb[4:-4] if pdb[:3]=="seq" else pdb[:-4]
print("WORKING on file: ",stem)
# make temporary directory 
try:
    subprocess.call("rm -r multi-pdb/; mkdir multi-pdb/",shell=True)
except:
    print("ERROR: generating multi-pdbs directory")
# pymol can not import multiple instances of the same pdb
# therefore create enough copies of target pdb to build brush
for n in np.arange(np.prod(size)):
    try:
        subprocess.call("cp {} multi-pdb/{}_{}.pdb".format(pdb,stem,n),
                        shell=True)
    except:
        print("ERROR: unable to copy to multi-pdb directory")

#########################################################
#               GENERATE PYMOL CMD FILE
######################################################### 

# load modules and pdbs into pymol
init_lines = "import pymol\nimport numpy as np\ncmd.delete('all')\n\
for i in np.arange({}):\n    cmd.load('multi-pdb/{}_{{}}.pdb'.format(i))\n\
cmd.orient('all')\n".format(np.prod(size), stem)
with open('temp.py','w') as f:
    f.write(init_lines)
# calculate the distance between "cylinders" of peptides
d = 2*args.rad+args.sep # / Angstroms
# create list of coordinates in the correct order to ensure sensible numbering
df = pd.DataFrame(list(product(np.arange(size[0]),np.arange(size[1]),
                                        repeat=1)))
df.sort_values(0, axis=0, ascending=False, inplace=True)
df.insert(value = np.zeros(np.prod(size)), column='Z', loc=0)
# pymol translate lines move imported pdbs to correct coordinates
for line in np.arange(np.prod(size)):
# multiply by d to get Angstrom values
    z,x,y = list(df.iloc[line].multiply(d))[:]
    translate_line = "cmd.translate([{}, {}, {}], \"{}_{}\" )\n"\
            .format(z,x,y,stem,line)
    with open('temp.py','a') as f:
        f.write(translate_line)
# pymol alter to renumber residues into continuous series
for i in np.arange(np.prod(size))[1:]:
    alter_line = "cmd.alter(\"{}_{}\" , 'resi=str(int(resi)+{l})')\n"\
            .format(stem,i,l=(args.len+1)*i if "_Ncapped" in args.pdb \
            else (args.len*i))
    with open('temp.py','a') as f:
        f.write(alter_line)
# save output pdb with detailed naming in new directory
outname = "brush_{}_{}x{}_{}A".format(stem,size[0],size[1],args.sep)
with open('temp.py','a') as f:
    f.write("cmd.extract('tosave', 'all')\ncmd.save('{o}/{o}.pdb', 'tosave')"\
            .format(o=outname))

#########################################################
#                       RUN PYMOL
######################################################### 

try:
    subprocess.call("mkdir {}".format(outname),shell=True)
    subprocess.call("pymol -cqi temp.py",shell=True)
except:
    print("ERROR: PyMol command failed.")

#########################################################
#               GENERATE GMX MAKE SCRIPT
######################################################### 

try:
    subprocess.call("cp -r amber_disp.ff {}/".format(outname),shell=True)
    subprocess.call("cp prep.mdp {}/".format(outname),shell=True)
except:
    print('ERROR:copy')


ter_search = "HA3 GLY * "
ter_command = "sed -i -e \"/{}$num /a TER\" ${{name}}.pdb".format(ter_search)
ri = [1,10,11,19]
ri_command ="ri {n[0]} | ri {n[1]}-{n[2]} | ri {n[3]}"\
        .format(n=[x+1 for x in ri] if "_Ncapped" in args.pdb else ri)

dim_xy = [ 2*N*(args.rad/10) + N*(args.sep/10) for N in size]
dim_z = np.ceil(0.35 * args.len + 2.4 + 1)
box_command = ('echo -e \"System\n System\" | $GMX editconf '
               '-f ${{name}}_b3.gro -o ${{name}}_box.gro '
               '-bt triclinic -box {} {} {} -n i.ndx'\
                .format(dim_xy[0],dim_xy[1],dim_z))

index_command = "echo -e \""
if "_Ncapped" in args.pdb:
    print("\n\nN Capping Detected\n\n")
    num_command = "num=$(($i * {} ))".format(args.len+1)
    for i in np.arange(np.prod(size)):
        index_command += " ri {} |".format((i+1)*(args.len+1))
else:
    num_command = "num=$(($i * {}))".format(args.len)
    for i in np.arange(np.prod(size)):
        index_command += " ri {} |".format((i+1)*args.len)

index_command += ("\\n 19 & 4\\n name 20 Anchor\\n q\" | "
                  "$GMX make_ndx -f ${name}.gro -o i.ndx")


with open('NEW_cubic.sh','r') as f:
    l= f.read()

l = l.replace("NAME",       "export name={}".format(outname))
l = l.replace("NUMPEP",     str(np.prod(size)))
l = l.replace("NUMCOMMAND", num_command)
l = l.replace("TERCOMMAND", ter_command)
l = l.replace("RICOMMAND", ri_command)
l = l.replace("BOXCOMMAND", box_command)
l = l.replace("NDXCOMMAND", index_command)

with open('{}/gmx_make.sh'.format(outname),'w') as f:
    f.write(l)


#########################################################
#               RUN GMX MAKE SCRIPT
######################################################### 

try:
    #subprocess.call("bash {}/gmx_make.sh".format(outname),shell=True)
    proc = subprocess.Popen("bash gmx_make.sh", cwd="{}".format(outname),
                            shell=True)
    proc.communicate()
except:
    print('ERROR:gmx')

