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


########################################################
#                   PARSING INPUTS
########################################################

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=d, epilog=" ")

# required arguments
parser.add_argument("-pdb", type=str, default='5ai0', help='PDB code (default: %(default)s)')
parser.add_argument("-fs", type=int, default=1, help='fs (default: %(default)s)')
parser.add_argument("-refpath", type=str, default='./reference.pdb', help='path to reference pdb for alignment in PyMol (default: %(default)s)')

# optional arguments
parser.add_argument("-gmx", type=str, default='/usr/local/gromacs/bin/gmx', help='Local Gromacs Executable Path (default: %(default)s)', required=False)
parser.add_argument("-pypath", type=str, default='/usr/bin/python', help='Local Python Executable Path (default: %(default)s)', required=False)
parser.add_argument("-compath", type=str, default='./IN_OUT_PDB/', help='Path for IN/OUT PDBs to make a comparison PyMOL session (default: %(default)s)', required=False)

# additional flags
#parser.add_argument("-prot", action="store_true", help='PROTonation of histidines (default: %(default)s)')

# define variable from user input arguments
args = parser.parse_args()
pdb     = args.pdb
fs      = args.fs
gmx_path = args.gmx
py_path = args.pypath
comp_path = args.compath

########################################################
#                  GENERATE OUT PDB
########################################################

# define file names and directory names 
wd = "{}_FS{}".format(pdb,fs)
stem = "{}-FS{}".format(pdb,fs)
colvar = "{}/{}.COLVAR.old".format(wd,stem)
xtc = "{}/{}_rw_fin.xtc".format(wd,stem)
tpr = "{}/{}_ext.tpr".format(wd,pdb)
out_name = "{}/{}_OUT.pdb".format(wd,stem)
# ensure that COLVAR file exists in the working directory
try:
    print("../UCB-350ns_Cterm/{}/{}.colvar ./{}".format(wd,pdb,colvar))
    subprocess.call('cp ../UCB-350ns_Cterm/{}/{}.colvar ./{}'.format(wd,pdb,colvar), shell=True)
except:
    print("ERROR: cannot find old COLVAR.")
    sys.exit()
# read in the old COLVAR file
with open(colvar) as f: 
    lines = f.readlines()
    data = [l for l in lines if l[0] not in ("@", "#")]
    data = [l.split()[:2] for l in data]
# extract time and pp.proj data
df = pd.DataFrame({'time':[float(x[0]) for x in data[:]],'proj':[float(x[1]) for x in data[:]]})
# sort by time so as to extract last suitable frame, i.e. greater difference from IN
df.sort_values('time',ascending=False, inplace=True)
# find the value where 3.4 < pp.proj < 3.5
df = df[ df['proj'].between(3.4,3.5,inclusive=True) ]
# get nearest timestamp to extract snapshot
timestamp = round(df['time'].iloc[0],-1)
# use Gromacs trjconv to extract the snapshot
try:
    subprocess.call("echo Protein_LIG | {} trjconv -s {} -f {} -o {} -b {t} -e {t} -n {}/i.ndx"\
            .format(gmx_path,tpr,xtc,out_name,wd,t=timestamp),shell=True)
except:
    print("ERROR: trjconv failed.")
    sys.exit()

########################################################
#               ALIGN ALL PDBS TO REF.
########################################################

# generate and run a simple PyMOL aligning script
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
    # write to a temporary file
    with open('temp.py','w') as f:
        f.write(init_lines)
        f.write("\n".join([cmd0,cmd1,cmd11,cmd2,cmd3]))
    # align using PyMOL
    try:
        subprocess.call("pymol -cqi temp.py",shell=True)
    except:
        print("ERROR: PyMol command failed.")
        sys.exit()
    # copy the pdbs to a comparison folder
    try:
        subprocess.call("cp {}_al.pdb {}".format(mobile_fn, comp_path), shell=True)
    except:
        print("ERROR: unable to copy to COMPARISON DIRECTORY")

########################################################
#           EDIT PDB ALIGN/ RMSD COLUMNS
########################################################

    # read in newly-aligned pdb file
    filename = "{}_al.pdb".format(mobile_fn)
    with open(filename) as f:
        lines = f.readlines()
        atom_data = [l for l in lines if l[0] not in ("@", "#")] 

    # reassign the final two columns ( ALIGN(Y/N) and RMSD(Y/N) )
    # for the ligand:               ALIGN = N, RMSD = Y
    # for backbone protein atoms:   ALIGN = Y, RMSD = N
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
    # save edited pdb
    ed_file = filename[:-4]+"_ed.pdb"
    with open(ed_file, 'w') as f: 
        for line in edited_data:
            f.write(line)

########################################################
#       CREATE PLUMED.DRIVER FILE & RUN DRIVER
########################################################

# generate plumed - Driver command file
with open("{}/plumed_4_driver.dat".format(wd), 'w') as f:
    f.write("rmsdI: RMSD REFERENCE={}/{}_IN_al_ed.pdb TYPE=OPTIMAL\n".format(wd,stem))
    f.write("rmsdO: RMSD REFERENCE={}/{}_OUT_al_ed.pdb TYPE=OPTIMAL\n".format(wd,stem))
    f.write("PRINT STRIDE=1 ARG=rmsdI.*,rmsdO.* FILE={}/{}.COLVAR.RW\n".format(wd,stem))
# call plumed Driver to generate new CV values
try:
    subprocess.call("plumed driver --mf_xtc {w}/{s}_rw_fin.xtc \
            --plumed {w}/plumed_4_driver.dat \
            --timestep 0.002 --trajectory-stride 5000"\
            .format(w=wd,s=stem),shell=True)
except:
    print("ERROR: Driver failed")
    sys.exit()


########################################################
#           COMBINE NEW AND OLD COLVAR
########################################################

# define COLVAR locations
old_COL_path="{}/{}.COLVAR.old".format(wd,stem)
new_COL_path="{}/{}.COLVAR.RW".format(wd,stem)
# read in old COLVAR header for column names
with open(old_COL_path) as f:
    head = f.readlines()[0]
head = head.split()[2:]
# read in old COLVAR file
old_COL = pd.read_csv(old_COL_path,delim_whitespace=True,names=head,skiprows=1)
# round timestamps to ensure successful merging
old_COL['int_time'] = old_COL['time'].astype(int)
# remove duplicate lines created by restarts
old_COL = old_COL.drop_duplicates(subset='int_time',keep='last')
# cutdown old COLVAR - select every 5th line
old_COL = old_COL.iloc[::5,:]
# read in new COLVAR header for column names
with open(new_COL_path) as f:
    head = f.readlines()[0]
head = head.split()[2:]
# read in new COLVAR file
new_COL = pd.read_csv(new_COL_path,delim_whitespace=True,names=head,skiprows=1)
# check if old/new COLVAR have been made/edited correctly
if old_COL.shape[0] != new_COL.shape[0]:
    print("ERROR: COLVARs of different dimenstions.")
    sys.exit()
# round timestamps to ensure successful merging
new_COL['int_time'] = new_COL['time'].astype(int)
# merge COLVARs on rounded timestamp & select necessary columns
comb_COL = pd.merge(new_COL,old_COL[['pp.proj','pp.ext','meta.bias','int_time']],on='int_time')
# reorder the columns for output 
column_order = ['time','pp.proj','pp.ext','meta.bias','rmsdI','rmsdO']
comb_COL = comb_COL[column_order]
# ensure a reasonable number of decimal places in output
comb_COL = comb_COL.round(8)
# define output COLVAR name
comb_COL_path = "{}/{}.COLVAR.combined".format(wd,stem)
# write new header for COLVAR
with open(comb_COL_path,'w') as f:
    f.write("#! FIELDS "+" ".join(list(comb_COL.columns.values))+"\n")
# write out new combined COLVAR
comb_COL.to_csv(comb_COL_path,sep=" ",header=False,index=False,mode='a')

########################################################
#                   RUN REWEIGHT.PY
########################################################

# define final output Free Energy name
outfile = "{}/{}.FES.RW".format(wd,stem)
# run reweighting python script
# may need to add some of these are user options
try:
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
    print("SUCCESS: {} Output Reweighted Free Energy - {}".format(outfile))
except:
    print("ERROR: reweight.py failed.")

