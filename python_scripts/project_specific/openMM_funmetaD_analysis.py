"""
===============================================================================
                    OPENMM PROJECT SPECIFIC ANALYSIS
===============================================================================

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from glob import glob
import scipy.stats
import subprocess
import MDAnalysis as mda
from MDAnalysis.analysis import rms

# Directory locations
STEM = "/media/rhys/Storage/jctc2_rhys_2022"
DATA_DIR = f"{STEM}/gpfs_data"
NWDA_DIR = f"{STEM}/New_data/funnel_fragment_paper"
FIGS_DIR = f"{STEM}/Figures"
DRBX_DIR = "/home/rhys/Dropbox/RESEARCH/AA_RHYS/BB_JCTC2"

# Experimental Values
EXP_VALS = {
    "3U5J": -7.65,
    "3U5L": -8.45,
    "4HBV": -6.33,
    "4LR6": -6.11,
    "4MEQ": -4.91,
    "4PCI": -6.99,
    "4UYD": -5.59,
    "2WI3": -4.71,
    "2XHT": -8.13,
    "2YK9": -6.95,
    "3OW6": -9.08,
    "3K99": -9.85,
    "2WI2": -4.71,
    "3EKO": -9.14,
    "1JVP": -7.91,
    "1JVP_B": -7.91,
    "1PXJ": -7.08,
    "1WCC": -4.71,
    "2VTA": -5.09,
    "2VTH": -5.35,
    "2VTM": -4.09,
    "3TIY": -6.51,
}

# PDB Ids for each of the 3 systems
SYSTEMS = {
    "BRD4": ["3U5L", "3U5J", "4PCI", "4HBV", "4LR6", "4UYD", "4MEQ"],
    "HSP90": ["3K99", "3EKO", "3OW6", "2XHT", "2YK9", "2WI2", "2WI3"],
    "CDK2": ["1JVP", "1JVP_B", "1PXJ", "3TIY", "2VTH", "2VTA", "1WCC", "2VTM"],
}

# Define the limits of each CV per system for processing bias files -> FES
#             sys -method        proj      ext
CV_LIMITS = {
    "BRD4-fun-metaD": [0.8, 3.2, 0.0, 0.80],
    "HSP90-fun-metaD": [0.0, 3.7, 0.0, 1.35],
    "CDK2-fun-metaD": [0.3, 4.2, 0.0, 1.10],
    "BRD4-fun-RMSD": [0.8, 3.2, 0.0, 1.25],
    "HSP90-fun-RMSD": [0.0, 3.7, 0.0, 1.50],
    "CDK2-fun-RMSD": [0.3, 4.2, 0.0, 2.50],
}

# Total time (ns) for the simulations
T_MAX = 2000

# Default file names for topology and trajectory
# -- SAME AS DOM --
filenames = {
    "fun-metaD": ["solute.pdb", "full_trj_solute.dcd"],
    "fun-RMSD": ["equilibrated.pdb", "trj.dcd"],
}

# Same rmsd threshold for all systems for bound
# -- SAME AS DOM --
bound_rmsd = 2.5

# Differing unbound states, due to length of funnel changing
# -- DOM = no unbound used --
unbound_rmsd = {"BRD4": 15, "HSP90": 20, "CDK2": 25}

# import G.O.A.T. --> _bound_check
# import G.O.A.T. --> _identify_recross


def find_rx():
    # define delay (ns) - i.e. no RX will be counted within ns from the start
    # -- SAME AS DOM (for fun-metaD) --
    delay = 0

    # create storage df
    data = pd.DataFrame(columns=["method", "system", "pdb", "rep", "number", "rx"])

    for method in ["fun-metaD", "fun-RMSD"]:
        if "RMSD" in method:
            location = f"{STEM}/Fun-RMSD/analysis/RMSD_Data"
            rep_range = [0, 1, 2]
        else:
            location = NWDA_DIR + "/RMSD_Data"
            rep_range = [3, 4, 5]
        for system in SYSTEMS.keys():
            for pdb in SYSTEMS[system]:
                for i in rep_range:
                    if method == "fun-metaD" and system == "BRD4" and i == 5:
                        continue
                    df = pd.DataFrame()
                    with open(f"{location}/{method}/{system}/{pdb}_{i}.p", "rb") as f:
                        rmsd_data = pickle.load(f)
                    df["RMSD"] = rmsd_data
                    # add time and convert to ns
                    # df['time'] = np.arange(0,len(df)) if 'RMSD' in method else np.arange(0,len(df))*0.1
                    df["time"] = np.arange(0, len(df))
                    # run recrossing counting function
                    N, rx = identify_recross(
                        df,
                        "RMSD",
                        bound=[bound_rmsd, 0.0],
                        unbound=unbound_rmsd[system],
                    )
                    # filter out initial values due to the defined delay
                    rx = [1.0] + [x for x in rx[1:] if x > delay]
                    # add values to data storage
                    data.loc[len(data.index)] = [method, system, pdb, i, N, rx]
    # save data
    data.to_hdf(f"{STEM}/NEW_rx_from_rmsd.h5", key="df", mode="w")


def _work_out_biasfile(bias_dir, t):
    """Calculate the bias file in a bias directory that comes after a certain
    time, t.
    (Entirely derived from Dom's method of finding the files, which uses
    the total bias in a file to calculate when it was made along the
    trajectory.)"""
    # Create a dictionary of the total bias in each bias file.
    d = {}
    for bias_file in glob(f"{bias_dir}/*"):
        bias = np.load(bias_file)
        total_bias = np.sum(bias)
        d[bias_file] = total_bias

    # Sort the files based on their total bias
    sort_d = sorted(d.items(), key=lambda x: x[1])
    sorted_bias_files = [i[0] for i in sort_d]
    # Count the total number of bias files
    N = len(sorted_bias_files)

    # Calculate the set of times corresponding (approx.) to when each bias
    #  file was made
    t_per_file = np.linspace(0, T_MAX, N)

    # Find the index of the file that comes after time t
    index_to_use = next(i for i, x in enumerate(t_per_file) if x >= t)
    # N.B.
    # bisect.bisect_left(t_per_file, t)  # is faster for a large list

    # Returns the full path of the bias file.
    return sorted_bias_files[index_to_use]


def _write_fes_from_bias(bias_path, system, method, out_path):
    """Write a FES.dat file from the .npy bias file.
    Based entirely on Dom's methodology.
    Requires CV bounds to be defined."""

    print(f"INFO: Writing FES for {bias_path}")
    # Define dT for bias calculation
    deltaT = 300 * (10 - 1)
    # Load the pickled bias file
    bias = np.load(bias_path, allow_pickle=True)
    # Define the CV space
    xticks = np.linspace(
        CV_LIMITS[f"{system}-{method}"][0],
        CV_LIMITS[f"{system}-{method}"][1],
        len(bias[0, :]),
    )
    yticks = np.linspace(
        CV_LIMITS[f"{system}-{method}"][2],
        CV_LIMITS[f"{system}-{method}"][3],
        len(bias[:, 0]),
    )
    # Write the FES file (same format as PLUMED)
    # fes[y,x]
    with open(out_path, "w") as f:
        # Header line
        f.write("#proj,cv2,bias(kj/mol)\n")
        for y in range(len(bias[:, 0])):
            for x in range(len(bias[0, :])):
                line = (
                    f"{xticks[x]:14.9f} "
                    f"{yticks[y]:14.9f} "
                    f"{-((300+deltaT)/deltaT)*bias[y,x]:14.9f}\n"
                )
                # f.write('%14.9f %14.9f %14.9f\n'% (xticks[x], yticks[y], -((300+deltaT)/deltaT)*bias[y,x]))
                f.write(line)
    print(f"INFO: Successfully written FES for {out_path}")


def extract_fes(mode="per_rx", stride=None):
    # Read in the rx information.
    rmsd_data = pd.read_hdf(
        "/home/rhys/Dropbox/RESEARCH/Project_JCTC2/NEW_rx_from_rmsd.h5", key="df"
    )

    for method in ["fun-RMSD"]:
        if "RMSD" in method:
            rep_range = [0, 1, 2]
        else:
            rep_range = [3, 4, 5]
        for system in SYSTEMS.keys():
            for pdb in SYSTEMS[system]:
                for i in rep_range:
                    if method == "fun-metaD" and system == "BRD4" and i == 5:
                        continue
                    wd = f"{NWDA_DIR}/{method}/{system}/{i}_output/{i}_output/{pdb}"

                    if mode == "per_rx":
                        # Get list of rx from rmsd data.
                        recrossings = rmsd_data.loc[
                            (rmsd_data.method == method)
                            & (rmsd_data.system == system)
                            & (rmsd_data.pdb == pdb)
                            & (rmsd_data.rep == i)
                        ].rx.values[0][1:]

                        print(f"Extracting FESs for {len(recrossings)} Rxs")
                        for i, rx in enumerate(recrossings):
                            bias_path = _work_out_biasfile(f"{wd}/bias_dir", rx)
                            _write_fes_from_bias(
                                bias_path, system, method, f"{wd}/rxFES_{i}_{rx}.dat"
                            )
                    elif stride:
                        print(f"Extracting FESs with STRIDE = {stride}")
                        dname = f"{wd}/tmp_fes/"
                        # Make the working directory
                        try:
                            subprocess.run(["mkdir", "-p", dname], check=True)
                        except subprocess.CalledProcessError as error:
                            print(
                                "Error code:",
                                error.returncode,
                                ". Output:",
                                error.output.decode("utf-8"),
                            )
                        t_stamps = np.linspace(0, 2000, stride)
                        for i, rx in enumerate(t_stamps):
                            bias_path = _work_out_biasfile(f"{wd}/bias_dir", rx)
                            _write_fes_from_bias(
                                bias_path,
                                system,
                                method,
                                f"{wd}/tmp_fes/rxFES_{i}_{int(rx)}.dat",
                            )

                    else:
                        print("Incorrect MODE or missing STRIDE")
                        break


#                            proj       ext
proj_basins = {
    "BRD4_B": [0.0, 17.0, 0.0, 5.0],  # A 1 = Narrow definition
    "BRD4_U": [25.0, 30.0, 0.0, 2.5],  # B
    "CDK2_B": [7.0, 15.0, 0.0, 2.5],  # A
    "CDK2_U": [30.0, 40.0, 0.0, 3.0],  # B
    "HSP90_B": [6.0, 14.0, 0.0, 3.0],  # A
    "HSP90_U": [30.0, 35.0, 0.0, 3.0],
}  # B

proj_basins_wide = {
    "BRD4_B": [8.0, 18.0, 0.0, 6.0],  # A 2 = Wide definitio
    "BRD4_U": [25.0, 30.0, 0.0, 3.0],  # B
    "CDK2_B": [5.0, 16.0, 0.0, 8.0],  # A
    "CDK2_U": [30.0, 40.0, 0.0, 3.0],  # B
    "HSP90_B": [8.0, 16.0, 0.0, 8.0],  # A
    "HSP90_U": [30.0, 35.0, 0.0, 3.0],
}  # B

#                          proj       rmsd
rmsd_basins = {
    "BRD4_B": [0.0, 17.0, 0.0, 5.0],  # A
    "BRD4_U": [25.0, 30.0, 0.0, 12.0],  # B
    "CDK2_B": [7.0, 15.0, 0.0, 5.0],  # A
    "CDK2_U": [30.0, 40.0, 0.0, 20.0],  # B
    "HSP90_B": [6.0, 14.0, 0.0, 5.0],  # A
    "HSP90_U": [30.0, 35.0, 0.0, 14.0],
}  # B

# volume corrections as calculated from Laco's notebook
vol_corr = {
    "BRD4": 2.901782970733069,
    "HSP90": 3.0003470683770175,
    "CDK2": 2.955066838544255,
}


def make_final_FES():
    for method in ["fun-metaD", "fun-RMSD"]:
        if "RMSD" in method:
            rep_range = [0, 1, 2]
        else:
            rep_range = [3, 4, 5]
        for system in SYSTEMS.keys():
            for pdb in SYSTEMS[system]:
                for i in rep_range:
                    if method == "fun-metaD" and system == "BRD4" and i == 5:
                        continue
                    wd = f"{NWDA_DIR}/{method}/{system}/{i}_output/{i}_output/{pdb}"
                    d = {}
                    for bias_file in glob(f"{wd}/bias_dir/*"):
                        bias = np.load(bias_file)
                        total_bias = np.sum(bias)
                        d[bias_file] = total_bias

                    sort_d = sorted(d.items(), key=lambda x: x[1])
                    sorted_bias_files = [i[0] for i in sort_d]
                    final_bias = sorted_bias_files[-1]
                    _write_fes_from_bias(
                        final_bias, system, method, f"{wd}/FES_2micro.dat"
                    )


def make_tables():
    # print dG values for all systems
    for method in ["fun-metaD", "fun-RMSD"]:
        if "RMSD" in method:
            rep_range = [0, 1, 2]
            basins = rmsd_basins
        else:
            rep_range = [3, 4, 5]
            basins = proj_basins

        to_csv = [
            "System,Ligand,Replica,NRx,Rx1,tRx1,Rx2,tRx2,Rx3,tRx3,RxMax,tRxMax,Exp\n"
        ]
        for system in SYSTEMS.keys():
            for pdb in SYSTEMS[system]:
                for i in rep_range:
                    new_line = f"{system},{pdb},{i},"
                    if method == "fun-metaD" and system == "BRD4" and i == 5:
                        continue
                    wd = f"{NWDA_DIR}/{method}/{system}/{i}_output/{i}_output/{pdb}"
                    nRx = len(glob(f"{wd}/rxFES_*"))
                    new_line += f"{nRx},"
                    nFES = min(len(glob(f"{wd}/rxFES_*")), 3)
                    for n in range(nFES):
                        fes_path = glob(f"{wd}/rxFES_{n}_*")[0]
                        print(fes_path)
                        dg = calculate_delta_g(
                            fes_path,
                            basins[f"{system}_B"],
                            basins[f"{system}_U"],
                            vol_corr[system],
                        )
                        new_line += f"{dg:.3f},{fes_path.split('/')[-1].split('_')[-1].split('.')[0]},"
                    new_line += "-," * 2 * max(0, (3 - nFES))
                    if nFES:
                        final_rx = glob(f"{wd}/rxFES_{nRx-1}_*")[0]
                        dg_final = calculate_delta_g(
                            final_rx,
                            basins[f"{system}_B"],
                            basins[f"{system}_U"],
                            vol_corr[system],
                        )
                        new_line += f"{dg_final:.3f},{final_rx.split('/')[-1].split('_')[-1].split('.')[0]},"
                    else:
                        new_line += "-,-,"
                    new_line += f"{EXP_VALS[pdb]}"
                    to_csv.append(new_line + "\n")

        with open(
            f"/home/rhys/Dropbox/RESEARCH/AA_RHYS/BB_JCTC2/Results_&_Figures/{method.lower()}_dg_values_NEWBASINS.csv",
            "w",
        ) as f:
            f.writelines(to_csv)


def make_average_tables():
    # print dG values for all systems
    for method in ["fun-metaD", "fun-RMSD"]:
        if "RMSD" in method:
            rep_range = [0, 1, 2]
            basins = rmsd_basins
        else:
            rep_range = [3, 4, 5]
            basins = proj_basins

        to_csv = [
            "System,Ligand,Total Number of Rx,3rd Rx dG Value,3rd RX Time,Final Rx dG Value,Final Rx Time,Experimental dG Value\n"
        ]
        for system in SYSTEMS.keys():
            for pdb in SYSTEMS[system]:
                new_line = f"{system},{pdb},"
                #
                stats = [[], [], [], [], []]
                for i in rep_range:
                    if method == "fun-metaD" and system == "BRD4" and i == 5:
                        continue
                    wd = f"{NWDA_DIR}/{method}/{system}/{i}_output/{i}_output/{pdb}"
                    print(wd)
                    nRx = len(glob(f"{wd}/rxFES_*"))
                    stats[0].append(nRx)
                    nFES = min(nRx, 3)
                    if nFES:
                        print(nFES)
                        fes_path = glob(f"{wd}/rxFES_{nFES-1}_*")[0]
                        dg = calculate_delta_g(
                            fes_path,
                            basins[f"{system}_B"],
                            basins[f"{system}_U"],
                            vol_corr[system],
                        )
                        stats[1].append(dg)
                        stats[2].append(
                            float(fes_path.split("/")[-1].split("_")[-1].split(".")[0])
                        )
                        final_rx = glob(f"{wd}/rxFES_{nRx-1}_*")[0]
                        dg_final = calculate_delta_g(
                            final_rx,
                            basins[f"{system}_B"],
                            basins[f"{system}_U"],
                            vol_corr[system],
                        )
                        stats[3].append(dg_final)
                        stats[4].append(
                            float(final_rx.split("/")[-1].split("_")[-1].split(".")[0])
                        )
                    else:
                        continue
                if any(stats[0]):
                    print(stats)
                    new_line += ",".join(
                        [
                            f"{np.asarray(s).mean():.2f} += {np.asarray(s).std():.2f}"
                            for s in stats
                        ]
                    )
                else:
                    new_line += "-," * 4
                new_line += f",{EXP_VALS[pdb]}"
                to_csv.append(new_line + "\n")

        with open(
            f"/home/rhys/Dropbox/RESEARCH/AA_RHYS/BB_JCTC2/Results_&_Figures/{method.lower()}_AVERAGE_values_NEWBASINS.csv",
            "w",
        ) as f:
            f.writelines(to_csv)


def make_paper_tables():
    # print dG values for all systems
    for method in ["fun-metaD", "fun-RMSD"]:
        if "RMSD" in method:
            rep_range = [0, 1, 2]
            basins = rmsd_basins
        else:
            rep_range = [3, 4, 5]
            basins = proj_basins

        to_csv = [
            "System,Ligand,Replica,No. Rx,dG Rx1,t Rx1,dG Rx2,t Rx2,dG Rx3,t Rx3,dG Final Rx,t Final Rx,dG Experimental\n"
        ]
        for system in SYSTEMS.keys():
            for pdb in SYSTEMS[system]:
                #
                stats = [[], [], [], []]
                for i in rep_range:
                    # Start a new line
                    new_line = f"{system},{pdb},{i},"
                    # Skip missing systems
                    if method == "fun-metaD" and system == "BRD4" and i == 5:
                        continue
                    wd = f"{NWDA_DIR}/{method}/{system}/{i}_output/{i}_output/{pdb}"
                    print(wd)
                    # Calculate total number of Rx
                    nRx = len(glob(f"{wd}/rxFES_*"))
                    # Add total Rx number to line
                    new_line += f"{nRx},"
                    # Only look at first 3 Rx for table
                    nFES = min(nRx, 3)
                    # Add data to line for first 3 Rx - dG and t
                    for n in range(nFES):
                        fes_path = glob(f"{wd}/rxFES_{n}_*")[0]
                        dg = calculate_delta_g(
                            fes_path,
                            basins[f"{system}_B"],
                            basins[f"{system}_U"],
                            vol_corr[system],
                        )
                        new_line += f"{dg:.3f},{fes_path.split('/')[-1].split('_')[-1].split('.')[0]},"
                        stats[n].append(dg)
                    # Fill line with dashes if < 3 Rx
                    new_line += "-," * 2 * max(0, (3 - nFES))
                    # Add the last Rx data to line - dG and t
                    if nFES:
                        final_rx = glob(f"{wd}/rxFES_{nRx-1}_*")[0]
                        dg_final = calculate_delta_g(
                            final_rx,
                            basins[f"{system}_B"],
                            basins[f"{system}_U"],
                            vol_corr[system],
                        )
                        new_line += f"{dg_final:.3f},{final_rx.split('/')[-1].split('_')[-1].split('.')[0]},"
                        stats[3].append(dg_final)
                    else:
                        new_line += "-,-,"
                    # Add Experimental value to end of line
                    new_line += f"{EXP_VALS[pdb]}"
                    to_csv.append(new_line + "\n")

                next_line = f"{system},{pdb},Mean,,"
                for s in stats:
                    if len(s) > 1:
                        next_line += (
                            f"{np.asarray(s).mean():.2f} += {np.asarray(s).std():.1f},,"
                        )
                    else:
                        next_line += "-,,"
                # Add Experimental value to end of line
                next_line += f"{EXP_VALS[pdb]}"
                to_csv.append(next_line + "\n")

        with open(
            f"/home/rhys/Dropbox/RESEARCH/AA_RHYS/BB_JCTC2/Results_&_Figures/{method.lower()}_PAPER_table.csv",
            "w",
        ) as f:
            f.writelines(to_csv)


def make_pd_hdf():
    # print dG values for all systems
    for method in ["fun-metaD", "fun-RMSD"]:
        if "RMSD" in method:
            rep_range = [0, 1, 2]
            basins = rmsd_basins
        else:
            rep_range = [3, 4, 5]
            basins = proj_basins

        k = 0
        for system in SYSTEMS.keys():
            for pdb in SYSTEMS[system]:
                for i in rep_range:
                    # Skip missing systems
                    if method == "fun-metaD" and system == "BRD4" and i == 5:
                        continue
                    wd = f"{NWDA_DIR}/{method}/{system}/{i}_output/{i}_output/{pdb}"
                    # print(wd)

                    # Calculate total number of Rx
                    nRx = len(glob(f"{wd}/rxFES_*"))
                    print(nRx)
                    N = []
                    dG = []
                    t = []
                    # Add data
                    for n in range(nRx):
                        fes_path = glob(f"{wd}/rxFES_{n}_*")[0]
                        dg = calculate_delta_g(
                            fes_path,
                            basins[f"{system}_B"],
                            basins[f"{system}_U"],
                            vol_corr[system],
                        )

                        N.append(n)
                        dG.append(dg)
                        t.append(
                            float(fes_path.split("/")[-1].split("_")[-1].split(".")[0])
                        )

                    # print(np.asarray([N, dG, t]).T)
                    # backbone:
                    inp = pd.DataFrame(
                        columns=["nrx", f"R{i}dGrx", f"R{i}trx"],
                        data=np.asarray([N, dG, t]).T,
                    ).set_index("nrx")
                    inp_l = pd.concat({pdb: inp}, axis=1)
                    inp_s = pd.concat({system: inp_l}, axis=1)

                    if k == 0:
                        print("First time --> Creating Files")
                        inp_s.to_hdf(f"{DRBX_DIR}/{method}_dGRX.h5", key="df")
                        k += 1
                        continue
                    # print('Further time --> Reading Files & Adding Data')
                    new = pd.read_hdf(f"{DRBX_DIR}/{method}_dGRX.h5", key="df")

                    if any([(mi == inp_s.columns)[0] for mi in new.columns]):
                        print("Updating values in DataFrame.")
                        new.update(inp_s)
                    else:
                        print("Adding new values to DataFrame.")
                        new = new.join(inp_s, how="outer")
                    # reorder columns
                    new = new.iloc[:, new.columns.sortlevel(0, sort_remaining=True)[1]]
                    new.to_hdf(f"{DRBX_DIR}/{method}_dGRX.h5", key="df")
                    k += 1


def make_dGoT_hdf(method):
    # print dG values for all systems
    if "RMSD" in method:
        rep_range = [0, 1, 2]
        basins = rmsd_basins
        tmp_folder = "tmp_fes"
        mode = "doms"
    else:
        rep_range = [3, 4, 5]
        basins = proj_basins
        tmp_folder = "sh_fes"
        mode = "plumed"

    k = 0
    for system in SYSTEMS.keys():
        print(system)
        for pdb in SYSTEMS[system]:
            for i in rep_range:
                # Skip missing systems
                if method == "fun-metaD" and system == "BRD4" and i == 5:
                    continue
                wd = f"{NWDA_DIR}/{method}/{system}/{i}_output/{i}_output/{pdb}/{tmp_folder}"
                dat_files = sorted(
                    glob(f"{wd}/*.dat"),
                    key=lambda x: int(x.split(".")[0].split("_")[-1]),
                )
                # print(dat_files)
                dG = []
                t = []
                for fes_path in dat_files:
                    # if mode == 'doms':
                    # fes_data = pd.read_table(fes_path, sep="\s+", header=0, names=['proj','ext','val'])
                    # else:
                    # fes_data = pd.read_table(fes_path, sep="\s+", names=['proj','ext','val', 'err1', 'err2'], comment='#')
                    dg = calculate_delta_g(
                        fes_path,
                        basins[f"{system}_B"],
                        basins[f"{system}_U"],
                        vol_corr[system],
                        mode=mode,
                    )
                    dG.append(dg)
                    t.append(fes_path.split(".")[0].split("_")[-1])

                inp = pd.DataFrame(
                    columns=["t", f"R{i}"], data=np.asarray([t, dG]).T
                ).set_index("t")
                inp_l = pd.concat({pdb: inp}, axis=1)
                inp_s = pd.concat({system: inp_l}, axis=1)

                if k == 0:
                    print("First time --> Creating Files")
                    inp_s.to_hdf(f"{DRBX_DIR}/{method}_dGoverT.h5", key="df")
                    k += 1
                    continue
                new = pd.read_hdf(f"{DRBX_DIR}/{method}_dGoverT.h5", key="df")

                if any([(mi == inp_s.columns)[0] for mi in new.columns]):
                    print("Updating values in DataFrame.")
                    new.update(inp_s)
                else:
                    print("Adding new values to DataFrame.")
                    new = new.join(inp_s, how="outer")
                # reorder columns
                new = new.iloc[:, new.columns.sortlevel(0, sort_remaining=True)[1]]
                new = new.astype(float)
                new.to_hdf(f"{DRBX_DIR}/{method}_dGoverT.h5", key="df")
                k += 1


def fes_per_rx():
    # plot FES per rx for all systems
    for method in ["fun-metaD", "fun-RMSD"]:
        if "RMSD" in method:
            rep_range = [0, 1, 2]
            basins = rmsd_basins
        else:
            rep_range = [3, 4, 5]
            basins = proj_basins
        for system in SYSTEMS.keys():
            for pdb in SYSTEMS[system]:
                for i in rep_range:
                    if method == "fun-metaD" and system == "BRD4" and i == 5:
                        continue
                    wd = f"{NWDA_DIR}/{method}/{system}/{i}_output/{i}_output/{pdb}"

                    nRx = len(glob(f"{wd}/rxFES_*"))

                    if nRx < 2:
                        continue
                    fig, ax = plt.subplots(1, nRx, figsize=(nRx * 8, 6))
                    plt.suptitle(f"FES with Rx: {system} + {pdb} Rep. {i}")
                    plt.subplots_adjust(top=0.85, right=0.915)
                    # nFES = min(len(glob(f"{wd}/rxFES_*")), 3)

                    # Establish the highest energy value for normalising cmap
                    max_vals = []
                    for n in range(nRx):
                        # Load in each FES in order of RX.
                        fes_path = glob(f"{wd}/rxFES_{n}_*")[0]
                        print(fes_path)
                        fes = np.loadtxt(fes_path)
                        # Only considering the energy, converted to kcal/mol
                        data = fes[:, 2] / 4.184
                        # Find the max, absolute, finite value
                        max_non_inf = np.amax(np.abs(data[np.isfinite(data)]))
                        max_vals.append(max_non_inf)
                        print("VMAX: ", max_non_inf)
                    # Identify the max max for the upper end of the cmap.
                    cmax = max(max_vals) + 1
                    print(f"using: {cmax}")

                    # Plot each FES subplot
                    i = 0
                    for n in range(nRx):
                        # Load in each FES in order of RX.
                        fes_path = glob(f"{wd}/rxFES_{n}_*")[0]
                        fes = np.loadtxt(fes_path)
                        # Convert distances to Angstroms
                        x = fes[:, 0] * 10
                        y = fes[:, 1] * 10
                        # Convert energy to kcal/mol
                        z = fes[:, 2] / 4.184
                        # Find the max finite value of the energy
                        max_non_inf = np.amax(z[np.isfinite(z)])
                        # Normalise cmaps across all subplots
                        z = z + (max(max_vals) - max_non_inf)
                        # Only needed if using graphics.fes
                        # z = z*4.184
                        # Get the time of RX from the file name
                        t = float(fes_path.split(".")[0].split("_")[-1])
                        # Define the set of contours
                        conts = np.arange(0.0, cmax + 1, 2.0)
                        # Interpolate the grid points for a smoother plot
                        N = 200
                        xi = np.linspace(x.min(), x.max(), N)
                        yi = np.linspace(y.min(), y.max(), N)
                        zi = scipy.interpolate.griddata(
                            (x, y), z, (xi[None, :], yi[:, None]), method="cubic"
                        )
                        # Plot the colour map
                        cmap = ax[i].contourf(
                            xi,
                            yi,
                            zi,
                            conts,
                            levels=15,
                            cmap="RdYlBu",
                            antialiased=True,
                        )
                        # Plot the contours over the colour map
                        ax[i].contour(
                            xi,
                            yi,
                            zi,
                            conts,
                            colors="k",
                            levels=15,
                            linewidths=0.5,
                            alpha=0.5,
                            antialiased=True,
                        )
                        # Plot formatting
                        ax[i].set_title(f"t = {t:.0f} ns")
                        ax[i].set_xlabel("CV1 / ")
                        ax[i].set_ylabel("CV2 / ")
                        # Add the basins as dashed boxes
                        if basins is not None:
                            b1 = plt.Rectangle(
                                (basins[f"{system}_B"][0], basins[f"{system}_B"][2]),
                                (basins[f"{system}_B"][1] - basins[f"{system}_B"][0]),
                                basins[f"{system}_B"][3],
                                ls="--",
                                fc="none",
                                ec="k",
                                lw=2.0,
                            )
                            ax[i].add_patch(b1)
                            b2 = plt.Rectangle(
                                (basins[f"{system}_U"][0], basins[f"{system}_U"][2]),
                                (basins[f"{system}_U"][1] - basins[f"{system}_U"][0]),
                                basins[f"{system}_U"][3],
                                ls="--",
                                fc="none",
                                ec="k",
                                lw=2.0,
                            )
                            ax[i].add_patch(b2)
                        i += 1

                    # Add the colour bar with label to the right of the plots
                    cax = plt.axes([0.98, 0.11, 0.01, 0.77])
                    cbar = plt.colorbar(
                        cmap,
                        cax=cax,
                        aspect=10,
                    )
                    # ticks=np.arange(0., cmax, 20))
                    cbar.set_label("Free Energy / kcal/mol", fontsize=10)
                    # Save the png
                    fig.savefig(
                        f"/home/rhys/Dropbox/RESEARCH/AA_RHYS/BB_JCTC2/Results_&_Figures/FES_per_RX_withBasins/{method}/{system}_{pdb}_{i}_FESperRX.png",
                        dpi=300,
                        bbox_inches="tight",
                    )


if __name__ == "__main__":
    # 1 - CALCULATE THE RX
    # calc_rmsd()
    # find_rx()

    # 2, 3, 4
    # extract_fes(mode='per_rx')

    # 5 -
    # make_tables()
    # make_average_tables()
    # make_paper_tables()
    # make_pd_hdf()

    # extract_fes(mode='fixed interval', stride=41)

    # make_dGoT_hdf('fun-metaD')
    make_dGoT_hdf("fun-RMSD")

    # make_final_FES()

    # fes_per_rx()
