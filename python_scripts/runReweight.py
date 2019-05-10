d="""
===============================================================================
                            REWEIGHTING SCRIPT
===============================================================================
"""
import subprocess
import argparse
import sys
import pymol
import pymol.cmd as cmd
import pandas as pd
import numpy as np

import glob

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
wd = "{}_FS{}".format(pdb,fs)
stem = "{}-FS{}".format(pdb,fs)

colvar = "{}/{}.COLVAR.old".format(wd,stem)
xtc = "{}/{}_rw_fin.xtc".format(wd,stem)
tpr = "{}/{}_ext.tpr".format(wd,pdb)

try:
    print("../UCB-350ns_Cterm/{}/{}.colvar ./{}".format(wd,pdb,colvar))
    subprocess.call('cp ../UCB-350ns_Cterm/{}/{}.colvar ./{}'.format(wd,pdb,colvar), shell=True)
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

for state in ["IN","OUT"]:
    # Tell PyMOL we don't want any GUI features.
    #pymol.pymol_argv = ['pymol', '-Qic']
    # Call the function below before using any PyMOL modules.
    #pymol.finish_launching()
    # load reference pdb for aligning
    init_lines = "import pymol\ncmd.delete('all')\n"
    cmd0 = "cmd.load(\"{}\", object=\"ref\")".format(args.refpath)
    # load the IN/OUT pdb
    mobile_fn = "{}/{}_{}".format(wd,stem,state)
    cmd1 = "cmd.load(\"{}.pdb\", object=\"mob\")".format(mobile_fn)
    cmd11="mobile = \"mob and resi 15-315\"\ntarget = \"ref and resi 15-315\""
    cmd2 = "cmd.align(mobile, target, object=\"aligned\")"
    cmd3 = "cmd.save(\"{}_al.pdb\", \"mob\")".format(mobile_fn)
    with open('temp.py','w') as f:
        f.write(init_lines)
        f.write("\n".join([cmd0,cmd1,cmd11,cmd2,cmd3]))

    try:
        subprocess.call("pymol -cqi temp.py",shell=True)
    except:
        print("ERROR: PyMol command failed.")

########################################################
#           EDIT PDB ALIGN/ RMSD COLUMNS
########################################################

    filename = "{}_al.pdb".format(mobile_fn)

    with open(filename) as f:
        lines = f.readlines()
        atom_data = [l for l in lines if l[0] not in ("@", "#")] 

    edited_data = []
    lig_res = atom_data[-2].split()[3]
    for line in atom_data:      
        if line.split()[0] in ["TER","END"]: pass
        elif line.split()[3] == lig_res:
            s = list(line) 
            s[-19] = "1"
            s[-25] = "0"
            line = "".join(s)
        elif line.split()[2] not in ["CA","C","N","O"]:
            s = list(line)            
            s[-25] = "0"
            line = "".join(s)
        edited_data.append(line)

    ed_file = filename[:-4]+"_ed.pdb"

    with open(ed_file, 'w') as f: 
        for line in edited_data:
            f.write(line)

########################################################
#       CREATE PLUMED.DRIVER FILE & RUN DRIVER
########################################################

with open("{}/plumed_4_driver.dat".format(wd), 'w') as f:
    f.write("rmsdI: RMSD REFERENCE={}/{}_IN_al_ed.pdb TYPE=OPTIMAL\n".format(wd,stem))
    f.write("rmsdO: RMSD REFERENCE={}/{}_OUT_al_ed.pdb TYPE=OPTIMAL\n".format(wd,stem))
    f.write("PRINT STRIDE=1 ARG=rmsdI.*,rmsdO.* FILE={}/{}.COLVAR.RW\n".format(wd,stem))

try:
    subprocess.call("plumed driver --mf_xtc {w}/{s}_rw_fin.xtc \
            --plumed {w}/plumed_4_driver.dat \
            --timestep 0.002 --trajectory-stride 5000"\
            .format(w=wd,s=stem),shell=True)
except:
    print("ERROR: Driver failed")


########################################################
#           COMBINE NEW AND OLD COLVAR
########################################################

old_COL_path="{}/{}.COLVAR.old".format(wd,stem)
new_COL_path="{}/{}.COLVAR.RW".format(wd,stem)

with open(old_COL_path) as f:
    head = f.readlines()[0]
head = head.split()[2:]

old_COL = pd.read_csv(old_COL_path,delim_whitespace=True,names=head,skiprows=1)
old_COL['int_time'] = old_COL['time'].astype(int)
old_COL = old_COL.drop_duplicates(subset='int_time',keep='last')
old_COL = old_COL.iloc[::5,:]

with open(new_COL_path) as f:
    head = f.readlines()[0]
head = head.split()[2:]

new_COL = pd.read_csv(new_COL_path,delim_whitespace=True,names=head,skiprows=1)

if old_COL.shape[0] != new_COL.shape[0]:
    print("ERROR: COLVARs of different dimenstions.")
    sys.exit()

new_COL['int_time'] = new_COL['time'].astype(int)

comb_COL = pd.merge(new_COL,old_COL[['pp.proj','pp.ext','meta.bias','int_time']],on='int_time')
column_order = ['time','pp.proj','pp.ext','meta.bias','rmsdI','rmsdO']
comb_COL = comb_COL[column_order]
comb_COL = comb_COL.round(6)

comb_COL_path = "{}/{}.COLVAR.combined".format(wd,stem)
with open(comb_COL_path,'w') as f:
    f.write("#! FIELDS "+" ".join(list(comb_COL.columns.values))+"\n")
comb_COL.to_csv(comb_COL_path,sep=" ",header=False,index=False,mode='a')

########################################################
#           COMBINE NEW AND OLD COLVAR
########################################################
py_path = "/usr/bin/python"
outfile = "{}/{}.FES.RW".format(wd,stem)
try:
#    subprocess.call("{} reweight.py -bsf 10 -fpref 5ai0_FS1/fes/fes_ -nf 35 -fcol 3 -colvar 5ai0_FS1/5ai0-FS1.COLVAR.combined -biascol 4 -rewcol 5 6 -v".format(py_path), shell=True)
    subprocess.call("{pp} reweight.py \
            -bsf 10 \
            -fpref {w}/fes/fes_ \
            -nf 35 \
            -fcol 3 \
            -colvar {cCp} \
            -biascol 4 \
            -rewcol 5 6 \
            -outfile {of} \
            -v".format(pp=py_path,w=wd,cCp=comb_COL_path,of=outfile),shell=True)
except:
    print("ERROR: reweight.py failed.")


