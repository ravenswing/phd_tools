"""
PLOTTING FOR HYDROPHOBIC BRUSH
"""

import mdtraj as md
import numpy as np
import matplotlib
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use("pdf")
# import os
# import glob
# import sys
# import seaborn as sns

###############################################################################
#                                   INPUTS
###############################################################################

# Source directory for the files on your local system
srcdir = "/home/rhys/Desktop/gpfs/rhys/Project/Peptides/Historic & Hydrophilic/G_hydrophobic_brush/data/"


# mol = 'brush_3HIS+GLYx3_4x4_0.0A'
# mol = 'brush_3HIS+GLYx3_Ncapped_4x4_0.0A'
figure_path = "/home/rhys/Desktop/gpfs/rhys/Project/Peptides/Historic & Hydrophilic/G_hydrophobic_brush/figures/Thesis/"
num_pep = 16
# default_x = np.linspace(0,2500000.0,num=25001)
plt.style.use("fivethirtyeight")

###############################################################################
#                               ROLAVGPLOT
###############################################################################


def rol_avgplot2(data, window_size, no_of_std, ylims, ylabel, figure_name, save):
    keys = list(data.keys())
    if num_pep == 9:
        g, axs = plt.subplots(3, 3, figsize=(30, 15))
        # reorder = [6,1,5,4,0,2,8,3,7]
        reorder = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    elif num_pep == 12:
        g, axs = plt.subplots(3, 4, figsize=(30, 15))
        reorder = [11, 7, 3, 4, 10, 1, 0, 2, 8, 5, 6, 9]
    axs = axs.ravel()
    # x = np.linspace(0,500000.0,num=25001)
    x = np.linspace(0, 1000000.0, num=10001)
    colours = [
        "xkcd:red",
        "xkcd:orange",
        "xkcd:vibrant purple",
        "xkcd:maroon",
        "xkcd:cerulean",
        "xkcd:deep magenta",
        "xkcd:teal",
        "xkcd:green",
        "xkcd:purple",
        "xkcd:grapefruit",
        "xkcd:forest green",
        "xkcd:indigo",
    ]
    shade = 0.4
    c = 0
    for i in range(len(keys)):
        nm = keys[i]
        for p in data[nm]:
            b = reorder[p]
            df = pd.DataFrame(data[nm][b])
            rolling_mean = df[:].rolling(window_size, center=True).mean()
            rolling_std = df[:].rolling(window_size, center=True).std()
            num = np.shape(x)[0]
            mean = np.reshape(rolling_mean.values, num)
            minim = np.reshape((rolling_mean + (rolling_std * no_of_std)).values, num)
            maxim = np.reshape((rolling_mean - (rolling_std * no_of_std)).values, num)
            axs[p].fill_between(
                [a / 1000 for a in x],
                minim,
                maxim,
                alpha=(shade / 3),
                facecolor="{}".format(colours[c]),
            )
            axs[p].plot(
                [a / 1000 for a in x],
                mean,
                alpha=shade,
                linewidth=2.0,
                color="{}".format(colours[c]),
                zorder=20,
            )
            axs[p].set_title("Peptide " + str(b + 1), fontsize=22)
            # axs[p].set_xlabel("Simulation Time / $\mu$ s")
            axs[p].set_xlabel("Simulation Time (ns)", fontsize=16)
            axs[p].set_ylabel(ylabel, fontsize=16)
            axs[p].set_ylim(ylims)
            axs[p].set_xlim(0.0, 500)
            axs[p].set_xlim(0.0, 1000)
            axs[p].tick_params(labelsize=12)
            # axs[p].legend(legend, loc=1)
            c += 1
        c = 0
        shade = 1
    plt.subplots_adjust(hspace=0.4)
    if save == True:
        g.savefig(figure_path + mol + "/" + figure_name + ".png")


##############################################################################
#                               PROCESSING HELICITY
##############################################################################


def processHelicity(srcdir, mol, frames):
    print("Processing Helicity")
    col = ["pep" + str(x + 1) for x in range(num_pep)]
    df = pd.DataFrame(columns=col)
    print(df.head())
    xtc = "{}{m}/05-MD/MD_final.xtc".format(srcdir, m=mol)
    topol = "{}{m}/05-MD/{m}_protein.gro".format(srcdir, m=mol)

    for chunk in md.iterload(xtc, 100, top=topol):
        in_data = np.empty((100, num_pep))
        for i in range(num_pep):
            peptide_chunk = chunk.atom_slice(
                chunk.top.select(
                    "resid {} to {}".format(i * pep_length, (i + 1) * pep_length - 1)
                )
            )
            dssp = md.compute_dssp(peptide_chunk, simplified=False)
            # helicity = (dssp == 'H') | (dssp == 'G') | (dssp == 'I')
            helicity = (dssp == "E") | (dssp == "B")
            in_data[:, i] = np.sum(helicity, axis=1) / 1.0
        df = df.append(pd.DataFrame(in_data, columns=col), ignore_index=True)
    df["tot_hel"] = df.sum(axis=1)
    print(df.head())
    print(df.shape)
    # df.iloc[:frames, :].to_csv(mol+'_hel_data.csv')
    df.iloc[:frames, :].to_csv(mol + "_strand_data.csv")


def processHBonds(mol, limit):
    print("Processing Hydrogen Bonds")
    arcdir = srcdir + "ARC_files/"
    # mode = 'salt_bridges'
    mode = "h_bonds"
    # Defines interval list from which file names are made
    segments = np.linspace(0, 950, 20)
    # loop over segmented data
    for i in range(len(segments)):
        lower, upper = int(segments[i]) + 1, int(segments[i]) + 50
        target = "{m}_{l}-{u}_{z}.npy".format(m=mol, l=lower, u=upper, z=mode)
        print("Processing: " + target)
        chunk = np.load("{s}{m}/{t}".format(s=arcdir, m=mol, t=target, z=mode))

        binary_contacts = chunk < limit
        df = pd.DataFrame(np.sum(binary_contacts, 1), columns=["tot_hyb"])
        if i == 0:
            df.to_csv(mol + "_hyb_data.csv", header=True, index=False)
        else:
            df.to_csv(mol + "_hyb_data.csv", mode="a", header=None, index=False)


def processMetric(srcdir, mol, frames, metric):
    print("Processing " + metric)
    col = ["pep" + str(x + 1) for x in range(num_pep)]
    df = pd.DataFrame(columns=col)
    print(df.head())
    xtc = "{}{m}/05-MD/MD_final.xtc".format(srcdir, m=mol)
    topol = "{}{m}/05-MD/{m}_protein.gro".format(srcdir, m=mol)

    for chunk in md.iterload(xtc, 100, top=topol):
        in_data = np.empty((100, num_pep))
        for i in range(num_pep):
            peptide_chunk = chunk.atom_slice(
                chunk.top.select(
                    "resid {} to {}".format(i * pep_length, (i + 1) * pep_length - 1)
                )
            )
            # dssp = md.compute_dssp(peptide_chunk, simplified=False)
            if metric == "RGYR":
                in_data[:, i] = md.compute_rg(peptide_chunk)
            elif metric == "RMSD":
                in_data[:, i] = md.rmsd(peptide_chunk, peptide_chunk[0])
            elif metric == "E2E":
                calphas = peptide_chunk.top.select_atom_indices(selection="alpha")
                pair = (calphas[0], calphas[-1])
                mpair = np.matrix(pair)
                in_data[:, i] = np.concatenate(
                    md.compute_distances(peptide_chunk, mpair)
                )

        df = df.append(pd.DataFrame(in_data, columns=col), ignore_index=True)
    df["tot_hel"] = df.mean(axis=1)
    print(df.head())
    print(df.shape)
    df.iloc[:frames, :].to_csv(mol + "_" + metric + "_data.csv")


###############################################################################
#                               PLOTTING & RUNNING
###############################################################################
def run(processing, mol):
    """run the analysis and/or plot"""
    if processing[0]:
        processHelicity(srcdir, mol, 50000)
    elif processing[1]:
        processHBonds(mol, 0.28)

    else:
        helicity_data = pd.read_csv(mol + "_hel_data.csv", usecols=["tot_hel"])
        helicity_data["tot_hel"] = (
            helicity_data["tot_hel"] / helicity_data["tot_hel"].max()
        ) * 100
        hbond_data = pd.read_csv(mol + "_hyb_data.csv", usecols=["tot_hyb"])
        # print(hbond_data.head())
        data = helicity_data.merge(hbond_data, left_index=True, right_index=True)
        # print(data.shape)

        fontpath = "/home/rhys/.fonts/iosevka/iosevka-term-regular.ttf"
        prop = font_manager.FontProperties(fname=fontpath)
        matplotlib.rcParams["font.family"] = prop.get_name()
        hyb_col = "#ffc000"
        hel_col = "#002C56"

        window_size = 500
        font_sizes = [40, 36]
        x = np.arange(50000)
        for d in ["hel", "hyb"]:
            data[d + "_rm"] = data["tot_" + d].rolling(window_size, center=True).mean()
            data[d + "_std"] = data["tot_" + d].rolling(window_size, center=True).std()
            data[d + "_min"] = data[d + "_rm"] + data[d + "_std"]
            data[d + "_max"] = data[d + "_rm"] - data[d + "_std"]
        data.to_csv(mol + "_data_out.csv")
        fig, ax1 = plt.subplots(figsize=(20, 16))

        ax1.plot(data["hel_rm"], color=hel_col, linewidth=4.0, zorder=21)
        ax1.fill_between(
            x, data["hel_min"], data["hel_max"], alpha=0.3, facecolor=hel_col, zorder=20
        )
        ax1.set_ylabel(
            "Percentage Helicity", fontweight="medium", fontsize=font_sizes[0]
        )
        ax1.set_ylim(-2.0, 100.0)
        ax1.set_xlabel("Time (ns)", fontweight="medium", fontsize=font_sizes[0])
        ax1.set_xlim(0.0, 50000.0)
        plt.yticks(fontweight="medium", fontsize=font_sizes[1])
        ax2 = ax1.twinx()

        ax2.plot(data["hyb_rm"], color=hyb_col, linewidth=4.0)
        ax2.fill_between(
            x, data["hyb_min"], data["hyb_max"], alpha=0.3, facecolor=hyb_col
        )
        ax2.set_ylabel(
            "No. Hydrogen Bonds", fontweight="medium", fontsize=font_sizes[0]
        )
        ax2.set_ylim(-0.1, 5.0)
        plt.yticks(fontweight="medium", fontsize=font_sizes[1])

        ns = np.linspace(0, 1000, 11, dtype="int")
        ts = np.linspace(0, 50000, 11)
        ax1.set_xticks(ticks=ts)
        ax1.set_xticklabels(labels=ns, fontweight="medium", fontsize=font_sizes[1])
        #    ax1.set_title('$(E_4K_4)_2$ Helical Content and Hydrogen Bond Formation', fontweight='medium',fontsize=font_sizes[0], pad=20)
        ax1.set_title(
            "$(EK)_5$ Helical Content and Hydrogen Bond Formation",
            fontweight="medium",
            fontsize=font_sizes[0],
            pad=20,
        )
        # ax1.set_title('AEAK... Helical Content and Hydrogen Bond Formation', fontweight='medium',fontsize=font_sizes[0], pad=20)

        fig.savefig(
            figure_path + mol + "_double_plot.png",
            bbox_inches="tight",
            transparent=True,
            dpi=300,
        )


def plot_brush_csv(csv, ylims, ylabel):
    data = pd.read_csv(csv)

    g, axs = plt.subplots(4, 4, figsize=(30, 25))
    axs = axs.ravel()
    x = np.linspace(0, 1000000.0, num=10001)
    colours = [
        "xkcd:red",
        "xkcd:orange",
        "xkcd:vibrant purple",
        "xkcd:maroon",
        "xkcd:cerulean",
        "xkcd:deep magenta",
        "xkcd:teal",
        "xkcd:green",
        "xkcd:purple",
        "xkcd:grapefruit",
        "xkcd:forest green",
        "xkcd:indigo",
        "xkcd:navy",
        "xkcd:gray",
        "xkcd:yellow",
        "xkcd:pink",
    ]
    shade = 0.4
    c = 0
    window_size = 200
    x = np.arange(100000)
    for i in range(16):
        pep = "pep" + str(i + 1)
        datum = pd.DataFrame(data[pep], columns=[pep])
        if "hel" or "strand" in csv:
            print("USING HELICITY ADJUSTMENT")
            datum[pep] = (datum[pep] / (pep_length - 2)) * 100
        print(datum.head())
        datum[pep + "_rm"] = datum[pep].rolling(window_size, center=True).mean()
        datum[pep + "_std"] = datum[pep].rolling(window_size, center=True).std()
        datum[pep + "_min"] = datum[pep + "_rm"] + datum[pep + "_std"]
        datum[pep + "_max"] = datum[pep + "_rm"] - datum[pep + "_std"]

        axs[i].plot(
            datum[pep + "_rm"], linewidth=2.0, color="{}".format(colours[c]), zorder=20
        )
        axs[i].fill_between(
            x,
            datum[pep + "_min"],
            datum[pep + "_max"],
            alpha=(shade / 3),
            facecolor="{}".format(colours[c]),
        )

        axs[i].set_title("Peptide " + str(i + 1), fontsize=22)
        axs[i].set_xlabel("Simulation Time (ns)", fontsize=16)
        axs[i].set_ylabel(ylabel, fontsize=16)
        if "hel" or "strand" in csv:
            axs[i].set_ylim(0.0, 100.0)
        else:
            axs[i].set_ylim(ylims)
        axs[i].set_xlim(0.0, 50000.0)
        ns = np.linspace(0, 1000, 11, dtype="int")
        ts = np.linspace(0, 50000, 11)
        axs[i].set_xticks(ticks=ts)
        axs[i].set_xticklabels(labels=ns, fontsize=12)
        # axs[p].legend(legend, loc=1)
        c += 1
    plt.subplots_adjust(hspace=0.4)
    g.savefig(
        figure_path + csv + ".png", bbox_inches="tight", transparent=True, dpi=300
    )


def plot_indcomb(name, pep_length, ylims, ylabels, window_size):
    """plot the other things"""
    font_sizes = [32, 24]
    x = np.arange(100000)

    """
    capped_name = name[0]+"Ncapped_"+name[1]
    capped_data = pd.read_csv(capped_name, usecols=['tot_hel'])
    if "hel" or "strand" in name[1]:
        capped_data['tot_hel'] = (capped_data['tot_hel'] / (num_pep*((pep_length+1)-2)))*100
    capped_data.rename(columns={'tot_hel': 'capd'}, inplace=True)
    """
    uncapped_name = name[0] + name[1]
    uncapped_data = pd.read_csv(uncapped_name, usecols=["tot_hel"])
    if "hel" or "strand" in name[1]:
        uncapped_data["tot_hel"] = (
            uncapped_data["tot_hel"] / (num_pep * (pep_length - 2))
        ) * 100
    uncapped_data.rename(columns={"tot_hel": "uncapd"}, inplace=True)

    # data = capped_data.merge(uncapped_data, left_index=True, right_index=True)
    data = uncapped_data

    i = 0
    # for d in ['capd', 'uncapd']:
    for d in ["uncapd"]:
        data[d + "_rm"] = data[d].rolling(window_size, center=True).mean()
        data[d + "_std"] = data[d].rolling(window_size, center=True).std()
        data[d + "_min"] = data[d + "_rm"] + data[d + "_std"]
        data[d + "_max"] = data[d + "_rm"] - data[d + "_std"]

        fig, ax1 = plt.subplots(figsize=(20, 8))

        """
        ax1.plot(data[d+'_rm'], color=f'C{0+i}', linewidth=4.0, zorder=21)
        ax1.plot(data[d], color=f'C{0+i}', alpha=.5, linewidth=2.0, zorder=21)
        """

        ax1.plot(data[d + "_rm"], color="C2", linewidth=4.0, zorder=21)
        ax1.plot(data[d], color="C2", alpha=0.5, linewidth=2.0, zorder=21)

        ax1.set_ylabel(ylabels, fontweight="medium", fontsize=font_sizes[0])
        ax1.set_ylim(ylims)
        ax1.set_xlabel("Time (ns)", fontweight="medium", fontsize=font_sizes[0])
        ax1.set_xlim(0.0, 50000.0)
        ns = np.linspace(0, 1000, 11, dtype="int")
        ts = np.linspace(0, 50000, 11)
        ax1.set_xticks(ticks=ts)
        ax1.set_xticklabels(labels=ns, fontweight="medium", fontsize=font_sizes[1])
        ax1.legend(["Rolling Mean", "Raw Data"], fontsize=font_sizes[1])
        plt.yticks(fontweight="medium", fontsize=font_sizes[1])
        fig.savefig(
            figure_path + uncapped_name + "_" + d + ".png",
            bbox_inches="tight",
            transparent=True,
            dpi=300,
        )
        i += 1

    fig2, ax2 = plt.subplots(figsize=(20, 16))

    """
    i=0
    for d in ['capd', 'uncapd']:
        ax2.plot(data[d+'_rm'], color=f'C{0+i}', linewidth=4.0, zorder=21)
        ax2.fill_between(x, data[d+'_min'], data[d+'_max'], alpha=.3, facecolor=f'C{0+i}', zorder=20)
        i += 1
    """

    for d in ["uncapd"]:
        ax2.plot(data[d + "_rm"], color="C2", linewidth=4.0, zorder=21)
        ax2.fill_between(
            x, data[d + "_min"], data[d + "_max"], alpha=0.3, facecolor="C2", zorder=20
        )

    ax2.set_ylabel(ylabels, fontweight="medium", fontsize=font_sizes[0])
    ax2.set_ylim(ylims)
    ax2.set_xlabel("Time (ns)", fontweight="medium", fontsize=font_sizes[0])
    ax2.set_xlim(0.0, 50000.0)
    ns = np.linspace(0, 1000, 11, dtype="int")
    ts = np.linspace(0, 50000, 11)
    ax2.set_xticks(ticks=ts)
    ax2.set_xticklabels(labels=ns, fontweight="medium", fontsize=font_sizes[1])
    ax2.legend(["N-terminal Capped", "Uncapped"], fontsize=font_sizes[1])
    plt.yticks(fontweight="medium", fontsize=font_sizes[1])
    fig2.savefig(
        figure_path + uncapped_name + "_COMBIND.png",
        bbox_inches="tight",
        transparent=True,
        dpi=300,
    )


PROCESSING = [False, False]
# run(PROCESSING, mol)

# for mol in ['brush_3HIS+GLYx3_4x4_0.0A', 'brush_3HIS+GLYx3_Ncapped_4x4_0.0A']:
# for mol in ['brush_3HIS+GLYx7_4x4_0.0A',]:
# processHelicity(srcdir, mol, 100000)
# processMetric(srcdir, mol, 100000, "E2E")
# processMetric(srcdir, mol, 100000, "RGYR")
# processMetric(srcdir, mol, 100000, "RMSD")
"""
for mol in ['brush_3HIS+GLYx3_', 'brush_3HIS+GLYx3_Ncapped_']:
    pep_length = 13 if "Ncapped" in mol else 12
    print("Using PEP_LENGTH:  "+str(pep_length))
    plot_brush_csv(mol+'_strand_data.csv', [-0.2, 5.0], '% Strand (DSSP)')
    plot_brush_csv(mol+'_hel_data.csv', [-0.2, 7.0], '% Helicity (DSSP)')
"""

"""
plot_indcomb(["brush_3HIS+GLYx3_", "4x4_0.0A_hel_data.csv"], 12, [-2.0, 100.0], '% Helicity (DSSP)', 1000)
plot_indcomb(["brush_3HIS+GLYx3_", "4x4_0.0A_strand_data.csv"], 12, [-2.0, 100.0], '% Strand (DSSP)', 1000)
plot_indcomb(["brush_3HIS+GLYx3_", "4x4_0.0A_RMSD_data.csv"], 12, [0.0, 0.5], 'RMSD from linear / nm', 1000)
plot_indcomb(["brush_3HIS+GLYx3_", "4x4_0.0A_RGYR_data.csv"], 12, [0.0, 2.0], '$R_{gyr}$ / nm', 1000)
plot_indcomb(["brush_3HIS+GLYx3_", "4x4_0.0A_E2E_data.csv"], 12, [0.0, 5.0], 'End-to-end Distance / nm', 1000)

"""
plot_indcomb(
    ["brush_3HIS+GLYx7_", "4x4_0.0A_hel_data.csv"],
    28,
    [-2.0, 100.0],
    "% Helicity (DSSP)",
    1000,
)
plot_indcomb(
    ["brush_3HIS+GLYx7_", "4x4_0.0A_strand_data.csv"],
    28,
    [-2.0, 100.0],
    "% Strand (DSSP)",
    1000,
)
plot_indcomb(
    ["brush_3HIS+GLYx7_", "4x4_0.0A_RMSD_data.csv"],
    28,
    [0.0, 0.5],
    "RMSD from linear / nm",
    1000,
)
plot_indcomb(
    ["brush_3HIS+GLYx7_", "4x4_0.0A_RGYR_data.csv"],
    28,
    [0.0, 2.0],
    "$R_{gyr}$ / nm",
    1000,
)
plot_indcomb(
    ["brush_3HIS+GLYx7_", "4x4_0.0A_E2E_data.csv"],
    28,
    [0.0, 5.0],
    "End-to-end Distance / nm",
    1000,
)
