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

# local scripts
#import graphics as gr

########################################################
#                   PARSING INPUTS
########################################################
# take in user arguments such that the script can be called from a bash script
# and still be controlled

# initialise the argument parser
PARSER = argparse.ArgumentParser(\
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Reweighting Scripts", epilog=" ")

# required arguments
PARSER.add_argument("-pdb", type=str, default='5ai0',
                    help='PDB code for the system (default: %(default)s)')
PARSER.add_argument("-fs", type=int, default=1,
                    help='Funnel Side (1 or 2) (default: %(default)s)')
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
PARSER.add_argument("-nooutpdb", action="store_true",
                    help='Supply a manual OUT pdb')

# define variable from user input argument via the argument parser
ARGS = PARSER.parse_args()
PDB = ARGS.pdb
FS = ARGS.fs
GMX_PATH = ARGS.gmx
PY_PATH = ARGS.pypath
COMP_PATH = ARGS.compath
WEIGHTS = ARGS.weights

# define file names and directory names
wd = "{}_FS{}".format(PDB, FS)              # target working directory
stem = "{}-FS{}".format(PDB, FS)            # root for all filenames
# define paths to all necessary input files
colvar = "{}/{}_OLD.colvar".format(wd, stem)# output COLVAR from simulation
xtc = "{}/{}_dry_center.xtc".format(wd, stem)# output xtc from post-processing
tpr = "{}/mdrun.tpr".format(wd)             # tpr file used in the simulation
ndx = "{}/index_{}.ndx".format(wd,PDB)      # index file for the current system
out_name = "{}/{}_OUT.pdb".format(wd, stem) # OUT pdb name

# if using -download:
# store of username@address for ssh and rsync, keys = remote server name
SERVERS = {'archer':   'rhys@login.archer.ac.uk',
           'archer2':  'rhyse@login.archer.ac.uk',
           'cscs':     'revans@ela.cscs.ch',
           'thomas':   'zcqsrev@thomas.rc.ucl.ac.uk',
           'mn':       'cnio96742@dt01.bsc.es',
          }

########################################################
#                   DOWNLOAD FILES
########################################################
# download the files from the remote server the MD was run on.
# THIS DOES NOT WORK AND HAS NOT BEEN USED IN THIS ANALYSIS
# ALL DOWNLOADS HAVE BEEN DONE MANUALLY INTO THE PRODUCTION/ DIRECTORY
def download_from_server(server, remote_pth):
    """ download from remote server """
    remote_svr = SERVERS[server]
    ## make new directory and download COLVAR & HILLS files
    try:
        subprocess.call("mkdir ./monitor/{}_{}/".format(PDB, FS),
                        shell=True)
        subprocess.call("scp {}:{}{f}/{m}/{f}/COLVAR \
                        monitor/{m}_{f}/{m}-{f}_OLD.colvar"\
                        .format(remote_svr, remote_pth, f=FS, m=PDB),
                        shell=True)
        print("scp {}:{}{f}/{m}/{f}/COLVAR monitor/{m}_{f}/{m}.colvar"\
                .format(remote_svr, remote_pth, f=FS, m=PDB))
        subprocess.call("scp {}:{}{f}/{m}/{f}/HILLS \
                        monitor/{m}_{f}/{m}-{f}.hills"\
                        .format(remote_svr, remote_pth, f=FS, m=PDB),
                        shell=True)
    except:
        print('ERROR: Unable to sync COLVAR/HILLS files \
                from remote server.')

########################################################
#                      SUM_HILLS
########################################################
# use plumed sum_hills to generate:
# 1. the initial FES using the pp.proj and pp.ext CVs
# 2. the 35 FES files (10ns each) for convergence and reweighting in the
#    fes/ directory
def run_sumhills():
    """ run plumed sum hills """

    ## create a customised sumhills.sh to run both the sum_hills commands
    try:
        subprocess.call("echo \'#!/bin/bash/\' > ./sumhills.sh",
                        shell=True)
        subprocess.call("echo \'cd {w}/\' >> ./sumhills.sh".format(w=wd),
                        shell=True)
        # create the "OLD" fes file, based on the original CVs
        subprocess.call('echo \"plumed sum_hills \
                        --hills {s}.hills \
                        --outfile {s}_OLD.fes \
                        --mintozero \" >> ./sumhills.sh'\
                        .format(s=stem),
                        shell=True)
        # create the 35 fes files in the fes/ directory
        subprocess.call('echo \"plumed sum_hills \
                        --hills {s}.hills \
                        --kt 2.5 \
                        --stride 5000 \
                        --outfile fes/fes_ \
                        --mintozero \" >> ./sumhills.sh'\
                        .format(s=stem),
                        shell=True)
        subprocess.call("echo \'cd ../\' >> ./sumhills.sh",
                        shell=True)
    except:
        print('ERROR: Unable to run SUM_HILLS 1')
    # make the fes/ directory and run the command file just created
    try:
        subprocess.call("mkdir {w}/fes/".format(w=wd), shell=True)
        os.system('bash ./sumhills.sh')
    except:
        print('ERROR: sumhills2')

########################################################
#                  GENERATE OUT PDB
########################################################
# Uses a scoring system based on pp.proj, pp.ext and the minimum distance from
# the ligand and the protein to determine a suitable frame for the OUT.pdb for
# use when calculating the RMSDs for reweighting
# NB: this did not work perfectly for all systems first try, and some had to be
#     made manually
def generate_outpdb():
    """ generate OUT pdb file based on scoring system """
    # read in the old COLVAR file
    with open(colvar) as f:
        lines = f.readlines()
        # remove comment lines and extraneous columns
        data = [l for l in lines if l[0] not in ("@", "#")]
        data = [l.split()[:3] for l in data]
    # extract time and pp.proj and pp.ext data
    df = pd.DataFrame({'time':[float(x[0]) for x in data[:]],\
            'proj':[float(x[1]) for x in data[:]],\
            'ext':[float(x[2]) for x in data[:]]})
    # round timestamps to ensure successful merging
    df['int_time'] = df['time'].astype(int)
    # remove duplicate lines created by restarts
    df = df.drop_duplicates(subset='int_time', keep='last')
    # cutdown COLVAR to match traj - select every 5th line
    df = df.iloc[::5, :]
    print(df.head(), df.shape)
    # define name for minimum distance data file
    dist_xvg = "{}/{}_mindist.xvg".format(wd, stem)
    # utilise Gromacs mindist to get minimum distance between protein and ligand
    try:
        subprocess.call("echo 1 13 | {} mindist -f {} -s {} -od {} -n {}"\
                .format(GMX_PATH, xtc, tpr, dist_xvg, ndx), shell=True)
    except:
        print("ERROR: gmx mindist failed")
        sys.exit()
    # read in mindist output
    with open(dist_xvg, 'r') as f:
        lines = f.readlines()
        md = [float(l.split()[1]) for l in lines if l[0] not in ("@", "#")]
    # add mindist data to pd.DataFrame
    #print(len(md))
    df['min_dist'] = md
    # filter the data to within basic thresholds on proj, ext and min. distance
    df = df[
            (df.proj > 3.1) & (df.proj < 4.0) \
            & (df.ext < 0.5)
            & (df.min_dist > 1.2)
           ]\
            if not ARGS.nocut else \
         df[
            (df.proj > 3.1) & (df.proj < 4.0) \
            | (df.ext < 0.5) \
            | (df.min_dist > 1.2)
           ]
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

    try:
        subprocess.call("echo Protein_LIG | {} trjconv -s {} -f {} -o {} \
                        -b {t} -e {t} -n {ind}"\
                        .format(GMX_PATH, tpr, xtc, out_name,
                                t=timestamp, ind=ndx), shell=True)
    except:
        print("ERROR: trjconv failed.")
        sys.exit()

########################################################
#               ALIGN ALL PDBS TO REF.
########################################################
# Use PyMOL to align the IN and newly-generated OUT pdbs, both are aligned
# to a reference pdb specified in refpath (5ai0-FS1 from CT_Aligned dir.)
def align_pdbs():
    """ generate and run a simple PyMOL aligning script """
    # operates the same for both IN and OUT, aligning both to the reference
    for state in ["IN", "OUT"]:
            # FOR LOADING PYMOL WITHOUT CREATING A SCRIPT...
            # ... DIDN'T WORK ON ARC FOR SOME REASON
            # Tell PyMOL we don't want any GUI features.
            #pymol.pymol_argv = ['pymol', '-Qic']
            # Call the function below before using any PyMOL modules.
            #pymol.finish_launching()
        # load reference pdb for aligning
        init_lines = "import pymol\ncmd.delete('all')\n"
        # load the reference pdb as chosen by refpath
        cmd0 = "cmd.load(\"{}\", object=\"ref\")".format(ARGS.refpath)
        # load the IN/OUT pdb
        mobile_fn = "{}/{}_{}".format(wd, stem, state)
        cmd1 = "cmd.load(\"{}.pdb\", object=\"mob\")".format(mobile_fn)
        # select only the core C-terminal residues
        cmd2 = "mobile = \"mob and resi 15-315\"\ntarget = \"ref and resi 15-315\""
        # run the alignment, creating "aligned" object
        cmd3 = "cmd.align(mobile, target, object=\"aligned\")"
        # save the aligned structure to a pdb file
        cmd4 = "cmd.save(\"{}_al.pdb\", \"mob\")".format(mobile_fn)
        # write to a temporary file for PyMOL to execute
        with open('temp.py', 'w') as f:
            f.write(init_lines)
            f.write("\n".join([cmd0, cmd1, cmd2, cmd3, cmd4]))
        # call PyMOL to run the commands above
        try:
            subprocess.call("pymol -cqi temp.py", shell=True)
        except subprocess.CalledProcessError as e:
            print("ERROR: PyMol command failed. {}, {}"
                  .format(e.output, e.returncode))
            sys.exit()
        # copy the newly aligned pdb to a directory for comparison
        try:
            subprocess.call("cp {}_al.pdb {}".format(mobile_fn, COMP_PATH),
                            shell=True)
        except:
            print("ERROR: unable to copy to COMPARISON DIRECTORY")


########################################################
#           EDIT PDB ALIGN/ RMSD COLUMNS
########################################################
# edit the ALIGN and RMSD columns of the IN and OUT pdbs
def edit_pdb(in_origin, out_origin):
    """
    in_origin/out_origin = either "VMD" or "GMX"
    reassign the final two columns ( ALIGN(Y/N) and RMSD(Y/N) )
    for the ligand:               ALIGN = N, RMSD = Y
    for backbone protein atoms:   ALIGN = Y, RMSD = N
    for all other protein atoms:  ALIGN = N, RMSD = N
    """
    # set the origins of the pdb as specified in input
    # this is necessary due to the ALIGN and RMSD columns having different
    # default values when produced by VMD and Gromacs
    origins = {"IN":in_origin, "OUT":out_origin}
    for state in ["IN", "OUT"]:
        # define filename root and working filename
        mobile_fn = "{}/{}_{}".format(wd, stem, state)
        filename = "{}_al.pdb".format(mobile_fn)
        # read in newly-aligned pdb file
        with open(filename) as f:
            lines = f.readlines()
            atom_data = [l for l in lines if l[0] not in ("@", "#")]
        edited_data = []
        #  assign the residue number of the ligand from the end of the pdb
        lig_res = atom_data[-2].split()[3]
        #print(lig_res)
        #print([x.split()[0] in ["TER","END"] for x in atom_data[-40:-30]])

        # scan through and edit all the lines in the source pdb
        for line in atom_data:
            # remove all TER/END lines that are not an atom definition
            if line.split()[0] != "ATOM":
                continue
            # assign the ligand to be the only RMSD candidate and not aligned
            elif line.split()[3] == lig_res:
                s = list(line)
                s[-19] = "1"
                s[-25] = "0"
                line = "".join(s)
            # set only protein backbone atoms for alignment
            elif origins[state] == "GMX" and line.split()[2] not in ["CA", "C", "N", "O"]:
                s = list(line)
                s[-25] = "0"
                line = "".join(s)
            # set only protein backbone atoms for alignment
            elif origins[state] == "VMD" and line.split()[2] in ["CA", "C", "N", "O"]:
                s = list(line)
                s[-25] = "1"
                line = "".join(s)
            # add all changed lines to new list
            edited_data.append(line)
        # save edited pdb
        ed_file = filename[:-4]+"_ed.pdb"
        with open(ed_file, 'w') as f:
            for line in edited_data:
                f.write(line)

########################################################
#       CREATE PLUMED.DRIVER FILE & RUN DRIVER
########################################################
# use plumed Driver to create the reweighted COLVAR file
def driver_rmsd():
    """ CREATE PLUMED.DRIVER FILE & RUN DRIVER """
    # path for the reweighted COLVAR
    new_col_path = "{}/{}_REW.colvar".format(wd, stem)
    # generate plumed - Driver command file
    with open("{}/plumed_4_driver.dat".format(wd), 'w') as f:
        # RMSD IN
        f.write("rmsdI: RMSD REFERENCE={}/{}_IN_al_ed.pdb TYPE=OPTIMAL\n"\
                .format(wd, stem))
        # RMSD OUT
        f.write("rmsdO: RMSD REFERENCE={}/{}_OUT_al_ed.pdb TYPE=OPTIMAL\n"\
                .format(wd, stem))
        # print output to reweighted COLVAR
        f.write("PRINT STRIDE=1 ARG=rmsdI.*,rmsdO.* FILE={}\n"\
                .format(new_col_path))
    # call plumed Driver to generate new COLVAR
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
# clean up old COLVAR file and combine with newly generated, reweighted COLVAR
# (removes comment lines and removes duplicates timestamps)
def combine_colvars():
    """ combine old and reweighted colvars """
    # define COLVAR locations
    old_col_path = "{}/{}_OLD.colvar".format(wd, stem)
    new_col_path = "{}/{}_REW.colvar".format(wd, stem)
    # read in old COLVAR header for column names
    with open(old_col_path) as f:
        head = f.readlines()[0]
    head = head.split()[2:]
    # read in old COLVAR file into dataFrame
    # filters out comment lines and splits columns via whitespace
    old_col = pd.concat([df[df.time != "#!"] \
                        for df in pd.read_csv(old_col_path,
                                              delim_whitespace=True,
                                              names=head,
                                              skiprows=1,
                                              chunksize=1000)])
    # round timestamps to ensure successful merging
    old_col['int_time'] = old_col['time'].astype(float).astype(int)
    # remove duplicate lines created by restarts
    old_col = old_col.drop_duplicates(subset='int_time', keep='last')
    # cutdown old COLVAR to match trajectories by selecting every 5th line
    old_col = old_col.iloc[::5, :]

    # add every 10th line (and the second line) for GISMO colvar = 3503 lines
    gis_col = old_col.iloc[:2, :]
    gis_col = gis_col.append(old_col.iloc[10::10, :], ignore_index=True)
    # define path for the original GISMO COLVAR file
    gismo_col_path = "{}/{}_OLD_GISMO.colvar".format(wd, stem)
    # add the header line to the new COLVAR
    with open(gismo_col_path, 'w') as f:
        f.write("#! FIELDS "+" ".join(list(gis_col.columns.values))+"\n")
    # save the cutdown original COLVAR
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
    # define path for the new combined COLVAR file
    comb_col_path = "{}/{}_COMBND.colvar".format(wd, stem)
    # add the header line to the new COLVAR
    with open(comb_col_path, 'w') as f:
        f.write("#! FIELDS "+" ".join(list(comb_col.columns.values))+"\n")
    # save the combined COLVAR
    comb_col.to_csv(comb_col_path, sep=" ", header=False, index=False, mode='a')

    # add every 10th line (and the second line) for GISMO colvar = 3503 lines
    out_col = comb_col.iloc[:2, :]
    out_col = out_col.append(comb_col.iloc[10::10, :], ignore_index=True)
    # define path for the reweighted GISMO COLVAR file
    gismo_col_path = "{}/{}_REW_GISMO.colvar".format(wd, stem)
    # add the header line to the new COLVAR
    with open(gismo_col_path, 'w') as f:
        f.write("#! FIELDS "+" ".join(list(out_col.columns.values))+"\n")
    # save the cutdown reweighted COLVAR
    out_col.to_csv(gismo_col_path, sep=" ", header=False, index=False, mode='a')

########################################################
#                   RUN REWEIGHT.PY
########################################################
# run the reweight.py script made by Giorgio
def run_ext_script():
    """ run Giorgio's reweight.py script """
    # set path to reweighted COLVAR
    comb_col_path = "{}/{}_COMBND.colvar".format(wd, stem)
    # define final output Free Energy name
    outfile = "{}/{}_REW.fes".format(wd, stem)
    # run reweighting python script
    # may need to add some of these as user options
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
# produce a cutdown version of the initial trajectory for use with GISMO
def cutdown_traj():
    """ cutdown the trajectories using Gromacs trjconv  ready for GISMO """
    # define the output GISMO xtc path
    gismo_xtc_path = "{}/{}_GISMO.xtc".format(wd, stem)
    # call gmx trjconv with -dt 100 to cut down the trajectory
    try:
        subprocess.call("echo Backbone Protein_LIG | {} trjconv -s {} -f {} \
                        -o {} -fit rot+trans -dt 100 -n {}"\
                        .format(GMX_PATH, tpr, xtc,
                                gismo_xtc_path, ndx),
                        shell=True)
    except:
        print("ERROR: trjconv (GISMO) failed.")
        sys.exit()

########################################################
#               1-D FES & CONVERGENCE
########################################################
# create convergence FES files and 1-D fes for each CV
def one_d_fes():
    """ run plumed sum hills """
    # create a customised sumhills.sh to run both the sum_hills commands
    try:
        subprocess.call("echo \'#!/bin/bash/\' > ./conv.sh",
                        shell=True)
        subprocess.call("echo \'cd {w}/\' >> ./conv.sh".format(w=wd),
                        shell=True)
    except:
        print('ERROR: Unable to run SUM_HILLS 1')

    for cv_name in ['proj', 'ext']:
        try:
            subprocess.call("mkdir conv_{cv}/".format(cv=cv_name), shell=True)
            subprocess.call('echo \"plumed sum_hills \
                            --hills {s}.hills \
                            --kt 2.5 \
                            --outfile {s}_{cv}.fes \
                            --idw pp.{cv} \
                            --mintozero \" >> ./sumhills.sh'\
                            .format(s=stem, cv=cv_name),
                            shell=True)
            subprocess.call('echo \"plumed sum_hills \
                            --hills {s}.hills \
                            --kt 2.5 \
                            --stride 5000 \
                            --outfile conv_{cv}/fes_ \
                            --idw pp.{cv} \
                            --mintozero \" >> ./sumhills.sh'\
                            .format(s=stem, cv=cv_name),
                            shell=True)
        except:
            print('ERROR: Unable to make convergence script/dirs')
    try:
        subprocess.call("echo \'cd ../\' >> ./sumhills.sh",
                        shell=True)
        os.system('bash ./sumhills.sh')
    except:
        print('ERROR: sumhills2')


    # set path to reweighted COLVAR
    comb_col_path = "{}/{}_COMBND.colvar".format(wd, stem)
    # define final output Free Energy name
    rmsds = {'rmsdIN':5, 'rmsdOUT':6}
    for rmsd in rmsds:
        outfile = "{}/{}_{}.fes".format(wd, stem, rmsd)
        # run reweighting python script individually for each CV
        try:
            subprocess.call("{pp} reweight.py \
                    -bsf 10 \
                    -fpref {w}/fes/fes_ \
                    -nf 35 \
                    -fcol 3 \
                    -colvar {cCp} \
                    -biascol 4 \
                    -rewcol {rwc} \
                    -outfile {of} \
                    -v".format(pp=PY_PATH, w=wd, cCp=comb_col_path,
                               rwc=rmsds[rmsd], of=outfile),
                            shell=True)
            print("SUCCESS: Output Reweighted Free Energy - {}".format(outfile))
        except:
            print("ERROR: reweight.py failed.")


########################################################
#                 DELTA G CALCULATIONS
########################################################

def calculate_dg():
    dat_files = ['fes.dat',]
    for fes_file in dat_files:

        fes_data = pd.concat([df[df.time != "#!"] \
                            for df in pd.read_csv(fes_file,
                                                  delim_whitespace=True,
                                                  names=['cv', 'fes_val'],
                                                  skiprows=1,
                                                  chunksize=1000)])

########################################################
#                   EXECUTE ALL CODE
########################################################

# only run in this manner if the code is called not imported
if __name__ == '__main__':
    # only download if the -download flag is given
    if ARGS.download:
        download_from_server("mn", "/home/cnio96/cnio96742/scratch/UCB_metaD/")

    #run_sumhills()

    # PLOTTING THAT DOES NOT WORK YET
    #gr.hills_plot(PDB, FS)
    #gr.diffusion_plots(PDB, FS, 2)
    #gr.two_cv_contour(PDB, FS, 30)

    if not ARGS.nooutpdb:
        generate_outpdb()
    #align_pdbs()
    #edit_pdb("VMD", "VMD")
    #driver_rmsd()
    #combine_colvars()
    #run_ext_script()
    #cutdown_traj()
    one_d_fes()
