import mdtraj as md
import numpy as np
import pandas as pd

# import seaborn as sns
import matplotlib
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt


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
        g.savefig(FIGDIR + mol + "/" + figure_name + ".png")


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
    x = np.arange(50000)
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
    g.savefig(FIGDIR + csv + ".png", bbox_inches="tight", transparent=True, dpi=300)


def plot_indcomb(name, pep_length, ylims, ylabels, window_size):
    """plot the other things"""

    col = {"capd": "#d08770", "uncapd": "#b48ead"}
    font_sizes = [32, 24]
    x = np.arange(100000)

    uncapped_name = name[0] + name[1]
    # capped_name = name[0]+"Ncapped_"+name[1]

    # capped_data = pd.read_csv(capped_name, usecols=['tot_hel'])
    uncapped_data = pd.read_csv(uncapped_name, usecols=["tot_hel"])
    if "hel" or "strand" in name[1]:
        #    capped_data['tot_hel'] = (capped_data['tot_hel'] / (num_pep*((pep_length+1)-2)))*100
        uncapped_data["tot_hel"] = (
            uncapped_data["tot_hel"] / (num_pep * (pep_length - 2))
        ) * 100
    # capped_data.rename(columns={'tot_hel': 'capd'}, inplace=True)
    uncapped_data.rename(columns={"tot_hel": "uncapd"}, inplace=True)

    # data = capped_data.merge(uncapped_data, left_index=True, right_index=True)
    data = uncapped_data

    # for d in ['capd', 'uncapd']:
    for d in ["uncapd"]:
        data[d + "_rm"] = data[d].rolling(window_size, center=True).mean()
        data[d + "_std"] = data[d].rolling(window_size, center=True).std()
        data[d + "_min"] = data[d + "_rm"] + data[d + "_std"]
        data[d + "_max"] = data[d + "_rm"] - data[d + "_std"]

        fig, ax1 = plt.subplots(figsize=(20, 8))
        ax1.plot(data[d + "_rm"], color=col[d], linewidth=4.0, zorder=21)
        ax1.plot(data[d], color=col[d], alpha=0.5, linewidth=2.0, zorder=21)
        # ax1.fill_between(x, data[d+'_min'], data[d+'_max'],
        # alpha=.3, facecolor=col[d], zorder=20)
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
            FIGDIR + uncapped_name + "_" + d + ".png",
            bbox_inches="tight",
            transparent=True,
            dpi=300,
        )
