"""
===============================================================================
                            REWEIGHTING SCRIPT
===============================================================================
"""
import argparse
import os
import subprocess
import sys
import pandas as pd
import pymol
import pymol.cmd as cmd


########################################################
#                   PARSING INPUTS
########################################################

PARSER = argparse.ArgumentParser(\
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Reweighting Scripts", epilog=" ")

# required arguments
PARSER.add_argument("-pdb", type=str, default='5ai0',
                    help='PDB code (default: %(default)s)')
PARSER.add_argument("-fs", type=int, default=1,
                    help='fs (default: %(default)s)')
PARSER.add_argument("-refpath", type=str, default='./reference.pdb',
                    help='path to reference pdb for alignment in PyMol' +
                    '(default: %(default)s)')

# optional arguments
PARSER.add_argument("-gmx", type=str, default='/usr/local/gromacs/bin/gmx',
                    help='Local Gromacs Executable Path(default: %(default)s)',
                    required=False)
PARSER.add_argument("-pypath", type=str, default='/usr/bin/python',
                    help='Local Python Executable Path (default: %(default)s)',
                    required=False)
PARSER.add_argument("-compath", type=str, default='./IN_OUT_PDB/',
                    help='Path for IN/OUT PDBs to make a comparison PyMOL' +
                    'session (default: %(default)s)', required=False)
PARSER.add_argument("-weights", type=float, nargs='+',
                    help='Weights for [ Proj, Ext, Min Dist ] ' +
                    '(default: %(default)s)', required=False)

# additional flags
PARSER.add_argument("-nocut", action="store_true",
                    help='Apply strict cutoff to data (default: %(default)s)')
PARSER.add_argument("-download", action="store_true",
                    help='Download data from MN (default: %(default)s)')

# define variable from user input arguments
ARGS = PARSER.parse_args()
PDB = ARGS.pdb
FS = ARGS.fs
GMX_PATH = ARGS.gmx
PY_PATH = ARGS.pypath
COMP_PATH = ARGS.compath
WEIGHTS = ARGS.weights

# define file names and directory names
wd = "{}_FS{}".format(PDB, FS)
stem = "{}-FS{}".format(PDB, FS)
colvar = "{}/{}.COLVAR.old".format(wd, stem)
xtc = "{}/{}_rw_fin.xtc".format(wd, stem)
tpr = "{}/{}_ext.tpr".format(wd, PDB)
out_name = "{}/{}_OUT.pdb".format(wd, stem)

SERVERS = {'archer':   'rhys@login.archer.ac.uk',
           'archer2':  'rhyse@login.archer.ac.uk',
           'cscs':     'revans@ela.cscs.ch',
           'thomas':   'zcqsrev@thomas.rc.ucl.ac.uk',
           'mn':       'cnio96742@dt01.bsc.es',
          }

########################################################
#                   DOWNLOAD FILES
########################################################
def download_from_server(server, remote_pth):
    """ download from remote server """
    remote_svr = SERVERS[server]
    ## make new directory and download COLVAR & HILLS files
    try:
        subprocess.call("mkdir ./monitor/{}_{}/".format(PDB, FS),
                        shell=True)
        subprocess.call("scp {}:{}{f}/{m}/{f}/COLVAR \
                        monitor/{m}_{f}/{m}.colvar"\
                        .format(remote_svr, remote_pth, f=FS, m=PDB),
                        shell=True)
        print("scp {}:{}{f}/{m}/{f}/COLVAR monitor/{m}_{f}/{m}.colvar"\
                .format(remote_svr, remote_pth, f=FS, m=PDB))
        subprocess.call("scp {}:{}{f}/{m}/{f}/HILLS \
                        monitor/{m}_{f}/{m}.hills"\
                        .format(remote_svr, remote_pth, f=FS, m=PDB),
                        shell=True)
    except:
        print('ERROR: Unable to sync COLVAR/HILLS files \
                from remote server.')

########################################################
#                      SUM_HILLS
########################################################
def run_sumhills():
    """ run plumed sum hills """

    ## run SUM_HILLS to generate free energy surface
    try:
        subprocess.call("echo \'#!/bin/bash/\' > ./sumhills.sh",
                        shell=True)
        subprocess.call('echo \"plumed sum_hills \
                        --hills monitor/{m}_{f}/{m}.hills \
                        --outfile monitor/{m}_{f}/{m}.fes \
                        --mintozero \" >> ./sumhills.sh'\
                        .format(f=FS, m=PDB),
                        shell=True)
        os.system('bash sumhills.sh')
    except:
        print('ERROR: Unable to run SUM_HILLS')

########################################################
#                  GENERATE OUT PDB
########################################################
def generate_outpdb():
    """ generate out pdb file based on scoring system """
    # ensure that COLVAR file exists in the working directory
    try:
        print("../UCB-350ns_Cterm/{}/{}.colvar ./{}".format(wd, PDB, colvar))
        subprocess.call('cp ../UCB-350ns_Cterm/{}/{}.colvar ./{}'\
                .format(wd, PDB, colvar), shell=True)
    except:
        print("ERROR: cannot find old COLVAR.")
        sys.exit()
    # read in the old COLVAR file
    with open(colvar) as f:
        lines = f.readlines()
        data = [l for l in lines if l[0] not in ("@", "#")]
        data = [l.split()[:3] for l in data]
    # extract time and pp.proj data
    df = pd.DataFrame({'time':[float(x[0]) for x in data[:]],\
            'proj':[float(x[1]) for x in data[:]],\
            'ext':[float(x[2]) for x in data[:]]})
    # round timestamps to ensure successful merging
    df['int_time'] = df['time'].astype(int)
    # remove duplicate lines created by restarts
    df = df.drop_duplicates(subset='int_time', keep='last')
    # cutdown COLVAR to match traj - select every 5th lin
    df = df.iloc[::5, :]
    print(df.head(), df.shape)
    # utilies Gromacs mindist to get minimum distance between protein and ligand
    dist_xvg = "{}/{}_mindist.xvg".format(wd, stem)
    try:
        subprocess.call("echo 1 13 | {} mindist -f {} -s {} -od {} -n {}/i.ndx"\
                .format(GMX_PATH, xtc, tpr, dist_xvg, wd), shell=True)
    except:
        print("ERROR: gmx mindist failed")
        sys.exit()
    # read in mindist output
    with open(dist_xvg, 'r') as f:
        lines = f.readlines()
        md = [float(l.split()[1]) for l in lines if l[0] not in ("@", "#")]
    # add mindist data to pd.DataFrame
    df['min_dist'] = md
    # filter the data to within basic thresholds on proj, ext and min. distance
    df = df[(df.proj > 3.1) & (df.proj < 4.0) \
                & (df.ext < 0.5) & (df.min_dist > 1.2)] if not ARGS.nocut \
            else df[(df.proj > 3.1) & (df.proj < 4.0) \
            | (df.ext < 0.5) | (df.min_dist > 1.2)]
    # find the maximum min. distance value
    max_min = df['min_dist'].max()
    # calculate the score for determining "best" OUT frame
    def calculate_score(frame):
        """ OUT pdb scoring function """
        # optimal value of pp.proj = 3.5
        proj_score = WEIGHTS[0] * abs(3.5 - frame['proj'])
        # optimal value of pp.ext = 0
        ext_score = WEIGHTS[1] * abs(0.0 - frame['ext'])
        # optimal value of min. distance = greatest possible
        d_score = WEIGHTS[2] * (max_min - frame['min_dist'])
        # score is sum of differences between opt. values and values per frame
        return proj_score + ext_score + d_score
    # calculate the score for each remaining frame
    df['score'] = df.apply(calculate_score, axis=1)
    # timestamp of frame with lowest score, i.e. " best" frame
    timestamp = df.loc[df['score'].idxmin()].int_time
    print(timestamp, type(timestamp))

    #### OLD METHOD OF RANKING FRAMES ####
    # sort by time so as to extract last suitable frame, i.e. greater
    #df.sort_values('time',ascending=False, inplace=True)
    # find the value where 3.4 < pp.proj < 3.5
    #df = df[ df['proj'].between(3.4,3.5,inclusive=True) ]
    # get nearest timestamp to extract snapshot
    #timestamp = round(df['time'].iloc[0],-1)
    # use Gromacs trjconv to extract the snapshot

    try:
        subprocess.call("echo Protein_LIG | {} trjconv -s {} -f {} -o {} \
                        -b {t} -e {t} -n {}/i.ndx"\
                .format(GMX_PATH, tpr, xtc, out_name, wd, t=timestamp), shell=True)
    except:
        print("ERROR: trjconv failed.")
        sys.exit()

########################################################
#               ALIGN ALL PDBS TO REF.
########################################################

def align_pdbs():
    """ generate and run a simple PyMOL aligning script """
    for state in ["IN", "OUT"]:
        # Tell PyMOL we don't want any GUI features.
        #pymol.pymol_argv = ['pymol', '-Qic']
        # Call the function below before using any PyMOL modules.
        #pymol.finish_launching()
        # load reference pdb for aligning
        init_lines = "import pymol\ncmd.delete('all')\n"
        cmd0 = "cmd.load(\"{}\", object=\"ref\")".format(ARGS.refpath)
        # load the IN/OUT pdb
        mobile_fn = "{}/{}_{}".format(wd, stem, state)
        cmd1 = "cmd.load(\"{}.pdb\", object=\"mob\")".format(mobile_fn)
        cmd2 = "mobile = \"mob and resi 15-315\"\ntarget = \"ref and resi 15-315\""
        cmd3 = "cmd.align(mobile, target, object=\"aligned\")"
        cmd4 = "cmd.save(\"{}_al.pdb\", \"mob\")".format(mobile_fn)
        # write to a temporary file
        with open('temp.py', 'w') as f:
            f.write(init_lines)
            f.write("\n".join([cmd0, cmd1, cmd2, cmd3, cmd4]))
        # align using PyMOL
        try:
            subprocess.call("pymol -cqi temp.py", shell=True)
        except subprocess.CalledProcessError as e:
            print("ERROR: PyMol command failed. {}, {}"
                  .format(e.output, e.returncode))
            sys.exit()
        # copy the pdbs to a comparison folder
        try:
            subprocess.call("cp {}_al.pdb {}".format(mobile_fn, COMP_PATH),
                            shell=True)
        except:
            print("ERROR: unable to copy to COMPARISON DIRECTORY")


########################################################
#           EDIT PDB ALIGN/ RMSD COLUMNS
########################################################

def edit_pdb():
    """
    reassign the final two columns ( ALIGN(Y/N) and RMSD(Y/N) )
    for the ligand:               ALIGN = N, RMSD = Y
    for backbone protein atoms:   ALIGN = Y, RMSD = N
    """
    for state in ["IN", "OUT"]:
        mobile_fn = "{}/{}_{}".format(wd, stem, state)
        # read in newly-aligned pdb file
        filename = "{}_al.pdb".format(mobile_fn)
        with open(filename) as f:
            lines = f.readlines()
            atom_data = [l for l in lines if l[0] not in ("@", "#")]

        edited_data = []
        lig_res = atom_data[-2].split()[3]
        for line in atom_data:
            if line.split()[0] in ["TER", "END"]:
                pass
            elif line.split()[3] == lig_res:
                s = list(line)
                s[-19] = "1"
                s[-25] = "0"
                line = "".join(s)
            elif line.split()[2] not in ["CA", "C", "N", "O"]:
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

def driver_rmsd():
    """ CREATE PLUMED.DRIVER FILE & RUN DRIVER """
    # generate plumed - Driver command file
    with open("{}/plumed_4_driver.dat".format(wd), 'w') as f:
        f.write("rmsdI: RMSD REFERENCE={}/{}_IN_al_ed.pdb TYPE=OPTIMAL\n"\
                .format(wd, stem))
        f.write("rmsdO: RMSD REFERENCE={}/{}_OUT_al_ed.pdb TYPE=OPTIMAL\n"\
                .format(wd, stem))
        f.write("PRINT STRIDE=1 ARG=rmsdI.*,rmsdO.* FILE={}/{}.COLVAR.RW\n"\
                .format(wd, stem))
    # call plumed Driver to generate new CV values
    try:
        subprocess.call("plumed driver --mf_xtc {x} \
                --plumed {w}/plumed_4_driver.dat \
                --timestep 0.002 --trajectory-stride 5000"\
                .format(w=wd, x=xtc), shell=True)
    except:
        print("ERROR: Driver failed")
        sys.exit()

########################################################
#           COMBINE NEW AND OLD COLVAR
########################################################

def combine_colvars():
    """ combine old and reweighted colvars """
    # define COLVAR locations
    old_col_path = "{}/{}.COLVAR.old".format(wd, stem)
    new_col_path = "{}/{}.COLVAR.RW".format(wd, stem)
    # read in old COLVAR header for column names
    with open(old_col_path) as f:
        head = f.readlines()[0]
    head = head.split()[2:]
    # read in old COLVAR file
    old_col = pd.read_csv(old_col_path, delim_whitespace=True, names=head,
                          skiprows=1)
    # round timestamps to ensure successful merging
    old_col['int_time'] = old_col['time'].astype(int)
    # remove duplicate lines created by restarts
    old_col = old_col.drop_duplicates(subset='int_time', keep='last')
    # cutdown old COLVAR - select every 5th line
    old_col = old_col.iloc[::5, :]

    # cutdown and save GISMO COLVAR file from original COLVAR(3501 lines)
    gis_col = old_col.iloc[::10, :]
    gismo_col_path = "{}/{}.old.GIScolvar".format(wd, stem)
    with open(gismo_col_path, 'w') as f:
        f.write("#! FIELDS "+" ".join(list(gis_col.columns.values))+"\n")
    gis_col.to_csv(gismo_col_path, sep=" ", header=False, index=False, mode='a')

    # read in new COLVAR header for column names
    with open(new_col_path) as f:
        head = f.readlines()[0]
    head = head.split()[2:]
    # read in new COLVAR file
    new_col = pd.read_csv(new_col_path, delim_whitespace=True, names=head,
                          skiprows=1)
    # check if old/new COLVAR have been made/edited correctly
    if old_col.shape[0] != new_col.shape[0]:
        print("ERROR: COLVARs of different dimenstions.")
        sys.exit()
    # round timestamps to ensure successful merging
    new_col['int_time'] = new_col['time'].astype(int)
    # merge COLVARs on rounded timestamp & select necessary columns
    comb_col = pd.merge(new_col,
                        old_col[['pp.proj', 'pp.ext', 'meta.bias', 'int_time']],
                        on='int_time')
    # reorder the columns for output
    column_order = ['time', 'pp.proj', 'pp.ext', 'meta.bias', 'rmsdI', 'rmsdO']
    comb_col = comb_col[column_order]
    # ensure a reasonable number of decimal places in output
    comb_col = comb_col.round(8)

    comb_col_path = "{}/{}.COLVAR.combined".format(wd, stem)
    with open(comb_col_path, 'w') as f:
        f.write("#! FIELDS "+" ".join(list(comb_col.columns.values))+"\n")
    comb_col.to_csv(comb_col_path, sep=" ", header=False, index=False, mode='a')

    # cutdown and save GISMO COLVAR file from reweighted COLVAR (3501 lines)
    comb_col = comb_col.iloc[::10, :]
    gismo_col_path = "{}/{}.RW.GIScolvar".format(wd, stem)
    with open(gismo_col_path, 'w') as f:
        f.write("#! FIELDS "+" ".join(list(comb_col.columns.values))+"\n")
    comb_col.to_csv(gismo_col_path, sep=" ", header=False, index=False, mode='a')

########################################################
#                   RUN REWEIGHT.PY
########################################################

def run_ext_script():
    """ run Giorgio's reweight.py script """
    # define final output Free Energy name
    outfile = "{}/{}.FES.RW".format(wd, stem)
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
                -v".format(pp=PY_PATH, w=wd, cCp=comb_col_path, of=outfile),
                        shell=True)
        print("SUCCESS: Output Reweighted Free Energy - {}".format(outfile))
    except:
        print("ERROR: reweight.py failed.")

########################################################
#               CUTDOWN GISMO TRAJ.
########################################################

def cutdown_traj():
    """ cutdown the trajectories using Gromacs trjconv  ready for GISMO """
    gismo_xtc_path = "{}/{}_GISMO.xtc".format(wd, stem)

    try:
        subprocess.call("echo Backbone Protein_LIG | {} trjconv -s {} -f {} \
                        -o {} -fit rot+trans -dt 100 -n {}/i.ndx"\
                        .format(GMX_PATH, tpr, xtc,
                                gismo_xtc_path, wd, t=timestamp),
                        shell=True)
    except:
        print("ERROR: trjconv (GISMO) failed.")
        sys.exit()

########################################################
#                   EXECUTE ALL CODE
########################################################

if __name__ == '__main__':
    if ARGS.download:
        download_from_server("mn", "/home/cnio96/cnio96742/scratch/UCB_metaD/" )
    run_sumhills()
    generate_outpdb()
    align_pdbs()
    edit_pdb()
    driver_rmsd()
    combine_colvars()
    run_ext_script()
    cutdown_traj()
