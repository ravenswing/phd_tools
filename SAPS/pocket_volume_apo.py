import subprocess
import MDAnalysis as mda
from MDAnalysis.analysis import align
import traj_tools as tt

systems = {'a2b1': ['A769', 'PF739', 'SC4', 'MT47', 'MK87'],
           'a2b2': ['A769', 'PF739', 'SC4', 'MT47', 'MK87']}

DATA_DIR = '/home/rhys/Storage/ampk_metad_all_data/apos'
#TMPL_DIR = f"{DATA_DIR}/pockets/cut_templates"


def aligned_pdb(wd, pdb_name, ref_path):
    u = tt._init_universe(f"{wd}/{pdb_name}")
    protein = u.select_atoms("protein or resname S2P")
    with mda.Writer(f'{wd}/tmp_prot.pdb', protein.n_atoms) as W:
        for ts in u.trajectory:
            W.write(protein)

    mobile = tt._init_universe(f'{wd}/tmp_prot.pdb')
    ref = tt._init_universe(ref_path)
    aligner = align.AlignTraj(mobile, ref, select='backbone',
                              filename=f'{wd}/aligned.pdb').run()


def aligned_dcd(wd, xtc_name, top, ref_path):
    # ADD NFRAMES ARGUMENT AN LINSPACE FOR FRAME ITERATION!
    u = tt._init_universe([f"{wd}/{top}", f"{wd}/{xtc_name}"])
    protein = u.select_atoms("protein or resname S2P")
    with mda.Writer(f'{wd}/tmp_prot.xtc', protein.n_atoms) as W:
        for ts in u.trajectory[::200]:
            W.write(protein)

    mobile = tt._init_universe([f'{wd}/aligned.pdb', f'{wd}/tmp_prot.xtc'])
    ref = tt._init_universe(ref_path)
    aligner = align.AlignTraj(mobile, ref, select='backbone',
                              filename=f'{wd}/aligned.dcd').run()


def pocket_select(wd, out_name):
    mpck_cmd = ("mdpocket --trajectory_file aligned.dcd "
                "--trajectory_format dcd "
                f"-f aligned.pdb -o {out_name} -n 3.0")
    try:
        subprocess.run(mpck_cmd, cwd=wd, shell=True, check=True)
    except subprocess.CalledProcessError as error:
        print('Error code:', error.returncode,
              '. Output:', error.output.decode("utf-8"))


def pocket_volume(wd, out_name, ref_path):

    try:
        subprocess.run(f'cp {ref_path} {wd}', shell=True, check=True)
    except subprocess.CalledProcessError as error:
        print('Error code:', error.returncode,
              '. Output:', error.output.decode("utf-8"))

    mpck_cmd = ("mdpocket "
                "--trajectory_file aligned.dcd "
                "--trajectory_format dcd "
                f"-f aligned.pdb "
                f"--selected_pocket {ref_path.split('/')[-1]} "
                f"-o {out_name} "
                "-n 3.0 -v 10000")
    try:
        subprocess.run(mpck_cmd, cwd=wd, shell=True, check=True)
    except subprocess.CalledProcessError as error:
        print('Error code:', error.returncode,
              '. Output:', error.output.decode("utf-8"))


for system in systems.keys():

    out_dir = DATA_DIR
    ref = f"/home/rhys/Storage/ampk_metad_all_data/super_ref/{system}.pdb"

    for rep in ['R1', 'R2', 'R3']:
        wd = f"{DATA_DIR}/{system}/{rep}"

        """
        aligned_pdb(wd, f"{system}.pdb", ref)
        aligned_dcd(wd, f"{system}_apo_{rep}.nc", f"{system}.pdb", ref)

        try:
            subprocess.run('rm tmp_*', cwd=wd, shell=True, check=True)
        except subprocess.CalledProcessError as error:
            print('Error code:', error.returncode,
                  '. Output:', error.output.decode("utf-8"))

        pocket_select(wd, f"{system}_apo_{rep}")

        """
        pocket_volume(wd, f"{system}_{rep}_vol",
                        f"{DATA_DIR}/{system}_ADaM.pdb")
