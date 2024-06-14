"""
===============================================================================
                                    PLOTTING RUN

===============================================================================
"""

import glob
import subprocess
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# local files
import graphics as gr
import load_data as load

colours = [
    "#31859C",  # FS1 & BS1
    "#FFC000",  # Tunnel
    "#7030A0",  # FS2 & BS2
]


PLOT_ALL = False
PDB = "5akk"
FS = 1
# FS = None

if PLOT_ALL:
    data_dirs = glob.glob("5a*_FS*")
elif FS is None:
    data_dirs = [PDB + "_FS1", PDB + "_FS2"]
elif FS <= 2:
    data_dirs = ["{}_FS{}".format(PDB, FS)]
else:
    print("ERROR: not sure what to plot")
    sys.exit()


def all_plots(dir_list):
    """plot ALL THE THINGS"""
    for directory in dir_list:
        # extract system info.
        pdb, fs = directory.split("_")
        # make directory for plots
        subprocess.call("mkdir Figures/{}/".format(directory), shell=True)
        # Hills plots
        hills_data = load.hills("{}/{}-{}.hills".format(directory, pdb, fs))
        gr.hills_plot(hills_data, pdb, fs, "Figures/{}".format(directory))
        # Diffusion plots
        diff_data = load.colvar(
            "{}/{}-{}_OLD.colvar".format(directory, pdb, fs), "as_numpy"
        )
        gr.diffusion_plots(diff_data, pdb, fs, 2, "Figures/{}".format(directory))
        # OLD FES plots
        fes_data, axis_labels = load.fes(
            "{}/{}-{}_OLD.fes".format(directory, pdb, fs), False
        )
        gr.two_cv_contour(
            fes_data,
            pdb,
            fs,
            axis_labels,
            30,
            "OLD_FES",
            "Figures/{}".format(directory),
        )
        # REW FES plots
        fes_data, axis_labels = load.fes(
            "{}/{}-{}_REW.fes".format(directory, pdb, fs), True
        )
        gr.two_cv_contour(
            fes_data,
            pdb,
            fs,
            ["RMSD-IN", "RMSD-OUT"],
            30,
            "REW_FES",
            "Figures/{}".format(directory),
        )


def fes_multiplot(pdb_list, cbar_max):
    fig = plt.figure(figsize=(14, 14))
    axes = fig.subplots(3, 2, sharex=True, sharey=True)
    n = 0
    for pdb in pdb_list:
        f = 0
        for fs in ["FS1", "FS2"]:
            directory = "{}_{}".format(pdb, fs)
            plot_num = int("".join(["2", str(len(pdb_list)), str(n)]))
            print(plot_num)
            fes_data, axis_labels = load.fes(
                "{}/{}-{}_OLD.fes".format(directory, pdb, fs), False
            )
            cmap = gr.two_cv_contour(
                fes_data,
                pdb,
                fs,
                axis_labels,
                cbar_max,
                "OLD_FES",
                "Figures/{}".format(directory),
                axes[n, f],
            )
            f += 1
        n += 1

    plt.subplots_adjust(left=0.15)
    # Y LABELS
    A = np.arange(len(pdb_list))
    for i in A:
        fig.text(
            0.1,
            0.24 + (A[-i - 1] * 0.26),
            axis_labels[1] + " / nm",
            va="center",
            rotation="vertical",
            fontsize=10,
        )
        fig.text(0.06, 0.23 + (A[-i - 1] * 0.26), pdb_list[i], ha="center", fontsize=14)
    # X LABELS
    for i in [0, 1]:
        fig.text(
            0.34 + (i * 0.39), 0.065, axis_labels[0] + " / nm", ha="center", fontsize=10
        )
        fig.text(0.34 + (i * 0.39), 0.9, "FS" + str(i + 1), ha="center", fontsize=14)
    fig.subplots_adjust(hspace=0.05, wspace=0.05)
    fig.subplots_adjust(right=0.915)
    cax = plt.axes([0.93, 0.11, 0.01, 0.77])
    cbar = plt.colorbar(
        cmap, cax=cax, aspect=10, ticks=np.arange(0.0, cbar_max + 1, 2.0)
    )
    cbar.set_label("Free Energy / kcal/mol", fontsize=10)
    fig.savefig(
        "Figures/FES_multi_" + "_".join(pdb_list) + ".png", dpi=300, bbox_inches="tight"
    )


def fes_highlight(pdb):
    """Plot both Raw and Reweighted data in once muliplot for one system"""
    fig = plt.figure(figsize=(10, 8))
    axes = fig.subplots(2, 2, sharex="col", sharey="col")
    cbar_max = 20
    fes_list = ["OLD", "REW"]

    n = 0
    for fes_type in fes_list:
        f = 0
        is_rew = True if fes_type == "REW" else False
        for fs in ["FS1", "FS2"]:
            directory = "{}_{}".format(pdb, fs)
            plot_num = int("".join(["2", str(len(fes_list)), str(n)]))
            print(plot_num)
            fes_data, axis_labels = load.fes(
                "{}/{}-{}_{}.fes".format(directory, pdb, fs, fes_type), is_rew
            )
            cmap = gr.two_cv_contour(
                fes_data,
                pdb,
                fs,
                axis_labels,
                cbar_max,
                "OLD_{}".format(fes_type),
                "Figures/{}".format(directory),
                axes[f, n],
            )
            f += 1
        n += 1
    plt.subplots_adjust(left=0.15)
    labels = {
        "OLD": ["Raw Data", "pp.proj / nm"],
        "REW": ["Reweighted", axis_labels[0] + " / nm"],
    }
    # Y LABELS
    for i in [0, 1]:
        fig.text(0.04, 0.3 + (i * 0.38), "FS" + str(i + 1), ha="center", fontsize=14)
        fig.text(
            0.08,
            0.31 + (i * 0.38),
            "pp.ext / nm",
            va="center",
            rotation="vertical",
            fontsize=10,
        )

    # X LABELS
    for i in np.arange(len(fes_list)):
        fig.text(
            0.33 + (i * 0.41), 0.9, labels[fes_list[i]][0], ha="center", fontsize=14
        )
        fig.text(
            0.33 + (i * 0.41), 0.065, labels[fes_list[i]][1], ha="center", fontsize=10
        )

    fig.subplots_adjust(hspace=0.05, wspace=0.15)
    fig.subplots_adjust(right=0.915)
    cax = plt.axes([0.93, 0.11, 0.01, 0.77])
    cbar = plt.colorbar(
        cmap, cax=cax, aspect=10, ticks=np.arange(0.0, cbar_max + 1, 2.0)
    )
    cbar.set_label("Free Energy / kcal/mol", fontsize=10)
    plt.show()
    fig.savefig("Figures/FES_highlight_" + pdb + ".png", dpi=300, bbox_inches="tight")


def conv_multiplot(pdb_list, mode, cv="proj"):
    """Plot 4x2 convergence plots (one CV) for paper"""
    if mode == "all":
        fig = plt.figure(figsize=(8, 21))
        axes = fig.subplots(4, 2, sharex=True, sharey=True)
        n = 0
        for pdb in pdb_list:
            f = 0
            for fs in ["FS1", "FS2"]:
                lines = (
                    [300, 350, 400, 450, 500]
                    if "kk" not in pdb
                    else [300, 350, 400, 450, 490]
                )
                gr.convergence("{}_{}/conv_{}".format(pdb, fs, cv), lines, axes[n, f])
                f += 1
            n += 1
        plt.subplots_adjust(left=0.1)
        A = np.arange(len(pdb_list))
        for i in A:
            fig.text(
                0.04, 0.2 + (A[-i - 1] * 0.2), pdb_list[i], ha="center", fontsize=14
            )
        fig.savefig(
            "Figures/CONV_{}_multi_{}.png".format(cv, "_".join(pdb_list)),
            dpi=300,
            bbox_inches="tight",
        )
        # fig.savefig('Figures/convergence_plot.png', dpi=300, bbox_inches='tight')

    elif mode == "cut":
        cv_list = ["proj", "ext"]
        fig = plt.figure(figsize=(10, 12.5))
        axes = fig.subplots(6, 4, sharex="col", sharey=True)
        n = 0
        for pdb in pdb_list:
            f = 0
            for cv in cv_list:
                for fs in ["FS1", "FS2"]:
                    lines = (
                        [300, 350, 400, 450, 500]
                        if "kk" not in pdb
                        else [300, 350, 400, 450, 490]
                    )
                    gr.convergence(
                        "{}_{}/conv_{}".format(pdb, fs, cv), lines, axes[n, f]
                    )
                    f += 1
            n += 1

        c2p = {"5aly": 0, "5akk": 2, "5alp": 1, "5ai0": 2, "5ai5": 1, "5alt": 2}
        fig.subplots_adjust(left=0.15, hspace=0.2, wspace=0.1)
        A = np.arange(len(pdb_list))
        # Y LABELS
        for i in A:
            fig.text(
                0.06,
                0.16 + (A[-i - 1] * 0.132),
                pdb_list[i],
                ha="center",
                fontsize=14,
                color=colours[c2p[pdb_list[i]]],
            )
            fig.text(
                0.1,
                0.16 + (A[i] * 0.132),
                "$\Delta$G / kcal/mol",
                va="center",
                rotation="vertical",
                fontsize=10,
            )
        # X LABELS
        for i in np.arange(len(cv_list)):
            fig.text(
                0.33 + (i * 0.38), 0.92, "pp." + cv_list[i], ha="center", fontsize=14
            )
            for j in [0, 1]:
                n = 2 * i + j
                fig.text(
                    0.24 + (n * 0.19),
                    0.065,
                    "pp." + cv_list[i] + " / nm",
                    ha="center",
                    fontsize=10,
                )
                fig.text(
                    0.24 + (n * 0.19), 0.89, "FS" + str(j + 1), ha="center", fontsize=14
                )

        fig.legend(
            [str(x) + " nm" for x in lines] + ["Reweight"],
            loc="lower center",
            ncol=6,
            frameon=False,
        )
        fig.savefig("Figures/CONV_paper.png", dpi=300, bbox_inches="tight")
    else:
        print("Mode Undefined")


def dif_multiplot(pdb_list):
    """Plot 4x3 diffusion plots (not used)"""
    cv_list = ["proj", "ext"]
    fig = plt.figure(figsize=(10, 12.5))
    axes = fig.subplots(3, 4, sharex=True, sharey="col")
    n = 0
    for pdb in pdb_list:
        f = 0
        for cv in cv_list:
            for fs in ["FS1", "FS2"]:
                data = load.colvar(
                    "{p}_{f}/{p}-{f}_OLD.colvar".format(p=pdb, f=fs), "as_pandas"
                )
                gr.diffusion(data, cv, axes[n, f])
                f += 1
        n += 1

    c2p = {"5aly": 0, "5akk": 2, "5alp": 1, "5ai0": 2, "5ai5": 1, "5alt": 2}
    fig.subplots_adjust(left=0.15, hspace=0.2, wspace=0.3)
    A = np.arange(len(pdb_list))
    for i in A:
        fig.text(
            0.06,
            0.225 + (A[-i - 1] * 0.275),
            pdb_list[i],
            ha="center",
            fontsize=14,
            color=colours[c2p[pdb_list[i]]],
        )
        fig.text(
            0.1,
            0.225 + (A[i] * 0.275),
            "$\Delta$G / kcal/mol",
            va="center",
            rotation="vertical",
            fontsize=10,
        )

    for i in np.arange(len(cv_list)):
        fig.text(0.33 + (i * 0.38), 0.92, "pp." + cv_list[i], ha="center", fontsize=14)
        for j in [0, 1]:
            n = 2 * i + j
            fig.text(
                0.24 + (n * 0.19),
                0.065,
                "Simulation Time / ns",
                ha="center",
                fontsize=10,
            )
            fig.text(
                0.24 + (n * 0.19), 0.89, "FS" + str(j + 1), ha="center", fontsize=14
            )

    fig.savefig("Figures/DIFF_paper.png", dpi=300, bbox_inches="tight")


def dif_multiplot2(pdb_list):
    """Plot 6x2 diffusion plots for paper"""
    cv_list = ["proj", "ext"]
    c2p = {"5aly": 0, "5akk": 2, "5alp": 1, "5ai0": 2, "5ai5": 1, "5alt": 2}

    fig = plt.figure(figsize=(10, 14))
    axes = fig.subplots(6, 2, sharex=True, sharey="row")

    n = 0
    for pdb in pdb_list:
        c = 0
        for cv in cv_list:
            f = 0
            for fs in ["FS1", "FS2"]:
                data = load.colvar(
                    "{p}_{f}/{p}-{f}_OLD.colvar".format(p=pdb, f=fs), "as_pandas"
                )
                gr.diffusion(data, cv, axes[n + c, f], colours[c2p[pdb]])
                f += 1
            c += 1
        n += 2
    fig.subplots_adjust(left=0.15, hspace=0.2, wspace=0.05)
    A = np.arange(len(pdb_list))
    for i in A:
        fig.text(
            0.06,
            0.225 + (A[-i - 1] * 0.265),
            pdb_list[i],
            ha="center",
            fontsize=14,
            color=colours[c2p[pdb_list[i]]],
        )
        fig.text(
            0.1,
            0.3 + (A[-i - 1] * 0.265),
            "pp." + cv_list[0] + " / nm",
            va="center",
            rotation="vertical",
            fontsize=10,
        )
        fig.text(
            0.1,
            0.17 + (A[-i - 1] * 0.265),
            "pp." + cv_list[1] + " / nm",
            va="center",
            rotation="vertical",
            fontsize=10,
        )

    for i in np.arange(len(cv_list)):
        fig.text(0.33 + (i * 0.38), 0.9, "FS" + str(i + 1), ha="center", fontsize=14)
        fig.text(
            0.33 + (i * 0.38), 0.065, "Simulation Time / ns", ha="center", fontsize=10
        )
    # fig.legend([str(x)+' nm' for x in lines] + ['Reweight'], loc='lower center', ncol=6, frameon=False)
    fig.savefig("Figures/DIFF_paper2.png", dpi=300, bbox_inches="tight")


###############################################################################

pdb_orders = [
    ["5am1", "5am3", "5alg"],
    ["5alp", "5alh", "5ai5"],
    # ['5am4', '5aly', '5akh'],
    ["5am0", "5alx", "5aia"],
    ["5alo", "5alt", "5akg"],
    ["5akk", "5ai0", "5ak6"],
]
for pdb_list in pdb_orders:
    fes_multiplot(pdb_list, cbar_max=22)

# fes_highlight('5akk')

pdb_16 = [
    ["5am1", "5alg", "5alp", "5alh"],
    ["5ai5", "5aly", "5akh", "5am0"],
    ["5alx", "5aia", "5alo", "5alt"],
    ["5akg", "5akk", "5ai0", "5ak6"],
]
# for CV in ['proj','ext']:
# for pdb_list in pdb:
# conv_multiplot(pdb_list, 'all', cv=CV)

pdb_cut = ["5aly", "5akk", "5alp", "5ai0", "5ai5", "5alt"]
# conv_multiplot(pdb_cut, 'cut')

# dif_multiplot2(['5aly','5alp','5alt'])

plt.show()
