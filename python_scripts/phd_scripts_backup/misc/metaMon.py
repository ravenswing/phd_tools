d = """
===============================================================================
                                    MONITOR

===============================================================================
"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import os
import subprocess
import glob
import argparse
import sys

DOWNLOAD = True
# DOWNLOAD=False

SUMHILLS = True
old = False

########################################################
#                   SERVER SSH INFO
########################################################
servers = {
    "archer": "rhys@login.archer.ac.uk",
    "archer2": "rhyse@login.archer.ac.uk",
    "cscs": "revans@ela.cscs.ch",
    "thomas": "zcqsrev@thomas.rc.ucl.ac.uk",
}

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter, description=d, epilog=" "
)
parser.add_argument(
    "-ncv", type=int, default=1, help="input number of CVs (default: %(default)s)"
)
parser.add_argument(
    "-svr", type=str, default="archer", help="Remote Server (default: %(default)s)"
)
parser.add_argument(
    "-svp",
    type=str,
    default="/work/e280/e280/rhys/",
    help="Remote Server Path ; where all system directories are stored(default: %(default)s)",
)
parser.add_argument(
    "-s",
    type=str,
    default="/work/e280/e280/rhys/",
    help="System Name (default: %(default)s)",
)
parser.add_argument(
    "-d", type=str, default="06-RunMD", help="Directory Name (default: %(default)s)"
)
parser.add_argument(
    "-vmax",
    type=float,
    default="-100.0",
    help="ColorMap Vmax (default: %(default)s)",
    required=False,
)


########################################################
#                   PARSING INPUTS
########################################################
args = parser.parse_args()

ncv = args.ncv
remote_svr = servers[args.svr]
remote_pth = args.svp
mols = {
    args.s: [],
}
folder = args.d


########################################################
#                   DOWNLOAD FILES
########################################################
if DOWNLOAD:
    for mol in mols:
        ## make new directory and download COLVAR & HILLS files
        try:
            subprocess.call("mkdir ./monitor/{}_{}/".format(mol, folder), shell=True)
            subprocess.call(
                "scp {}:{}{m}/{f}/COLVAR monitor/{m}_{f}/{m}.colvar".format(
                    remote_svr, remote_pth, f=folder, m=mol
                ),
                shell=True,
            )
            subprocess.call(
                "scp {}:{}{m}/{f}/HILLS monitor/{m}_{f}/{m}.hills".format(
                    remote_svr, remote_pth, f=folder, m=mol
                ),
                shell=True,
            )
        except:
            print("ERROR: Unable to sync COLVAR/HILLS files from remote server.")
if SUMHILLS:
    for mol in mols:
        ## run SUM_HILLS to generate free energy surface
        try:
            subprocess.call("echo '#!/bin/bash/' > ./sumhills.sh", shell=True)
            subprocess.call(
                'echo "plumed sum_hills --hills monitor/{m}_{f}/{m}.hills --outfile monitor/{m}_{f}/{m}.fes --mintozero " >> ./sumhills.sh'.format(
                    f=folder, m=mol
                ),
                shell=True,
            )
            os.system("bash sumhills.sh")
        except:
            print("ERROR: Unable to run SUM_HILLS")


########################################################
#                   ANALYSE & PLOT
########################################################

plt.style.use("ggplot")
for mol in mols:
    ########################################################
    #                       FOR 1 CV
    ########################################################
    if ncv == 1:
        HILLS = mols.copy()
        HILLS[mol] = [[], []]
        filename = "./monitor/{m}_{f}/{m}.hills".format(f=folder, m=mol)
        with open(filename) as f:
            lines = f.readlines()
            data = [l for l in lines if l[0] not in ("@", "#")]
            data = [[float(val) for val in line.split()] for line in data]
            HILLS[mol][0] = [l[0] for l in data]
            HILLS[mol][1] = [l[3] for l in data]
        t = int(5 * round((HILLS[mol][0][-1] / 1000) / 5))
        plt.figure()
        plt.plot([x / 1000 for x in HILLS[mol][0]], HILLS[mol][1], label=mol)
        plt.legend()
        plt.title("{m}  |  {}  |  Heights  |  {}ns".format(folder, t, m=mol))
        plt.savefig(
            "monitor/{m}_{f}/{m}-{f}-Heights_{}.png".format(t, f=folder, m=mol),
            bbox_inches="tight",
            dpi=300,
        )

        FES = mols.copy()
        FES[mol] = [[], []]
        filename = "./monitor/{m}_{f}/{m}.fes".format(f=folder, m=mol)
        with open(filename) as f:
            lines = f.readlines()
            data = [l for l in lines if l[0] not in ("@", "#")]
            data = [[float(val) for val in line.split()[:2]] for line in data]
            FES[mol][0] = [l[0] for l in data]
            FES[mol][1] = [l[1] for l in data]
        plt.figure()
        plt.plot(FES[mol][0], [x / 4.184 for x in FES[mol][1]], label=mol)
        # plt.legend()
        plt.ylabel("FE (kcal/mol)")
        plt.xlabel(r"$\alpha$RMSD")
        plt.ylim([0, 20])
        plt.title("{m}".format(m=mol))
        plt.savefig(
            "monitor/{m}_{f}/{m}-{f}-FES_{}.png".format(t, f=folder, m=mol),
            bbox_inches="tight",
            dpi=300,
        )

        COL = mols.copy()
        filename = "./monitor/{m}_{f}/{m}.colvar".format(f=folder, m=mol)
        with open(filename) as f:
            lines = f.readlines()
            var = lines[0].split()[3]
            data = [l for l in lines if l[0] not in ("@", "#")]
            data = [[float(val) for val in line.split()] for line in data]
        COL[mol] = [[]] * len(data[0])
        for i in np.arange(len(data[0])):
            COL[mol][i] = [l[i] for l in data]
        plt.figure()
        plt.plot([x / 1000 for x in COL[mol][0]], COL[mol][1], label=mol)
        plt.legend()
        plt.title("{m}  |  {}  |  {}  |  {}ns".format(folder, var, t, m=mol))
        plt.savefig(
            "monitor/{m}_{f}/{m}-{f}_{}_{}.png".format(var, t, f=folder, m=mol),
            bbox_inches="tight",
            dpi=300,
        )

    ########################################################
    #                       FOR 2 CVs
    ########################################################
    if ncv == 2:
        vmax = args.vmax

        HILLS = mols.copy()
        HILLS[mol] = [[], []]
        filename = "./monitor/{m}_{f}/{m}.hills".format(f=folder, m=mol)
        with open(filename) as f:
            lines = f.readlines()
            data = [l for l in lines if l[0] not in ("@", "#")]
            data = [[float(val) for val in line.split()] for line in data]
            HILLS[mol][0] = [l[0] for l in data]
            HILLS[mol][1] = [l[5] for l in data]
        t = int(5 * round((HILLS[mol][0][-1] / 1000) / 5))
        plt.figure()
        plt.plot([x / 1000 for x in HILLS[mol][0]], HILLS[mol][1], label=mol)
        plt.legend()
        plt.title("{m}  |  {}  |  Heights  |  {}ns".format(folder, t, m=mol))
        plt.savefig(
            "monitor/{m}_{f}/{m}-{f}_Heights_{}.png".format(t, f=folder, m=mol),
            bbox_inches="tight",
            dpi=300,
        )

        ## read in and remove blanks from .fes file
        FES = mols.copy()
        FES[mol] = [[], [], []]
        filename = "./monitor/{m}_{f}/{m}.fes".format(f=folder, m=mol)
        with open(filename) as f:
            lines = f.readlines()
            data = [l for l in lines if l[0] not in ("@", "#")]
        data = [[float(val) for val in line.split()] for line in data]
        breaks = [i for i, e in enumerate(data) if e == []]
        # remove blank lines
        for index in sorted(breaks, reverse=True):
            del data[index]
        ## get number of bins and CV names from header
        n = int([l.split()[-1] for l in lines if "nbins_pp.proj" in l][0])
        [x_name, y_name] = lines[0].split()[2:4]
        ## organise data into plotable arrays
        split_data = [data[i : i + n] for i in range(0, len(data), n)]
        z = []
        for block in split_data:
            z.append([l[2] for l in block])
            FES[mol][1].append(block[0][1])
        FES[mol][0] = [l[0] for l in split_data[0]]
        FES[mol][2] = np.asarray(z)
        with open("./monitor/vmax.dat", "a+") as f:
            f.write("{} {} {}\r\n".format(mol, folder, np.amax(FES[mol][2])))

        x, y = np.meshgrid(FES[mol][0], FES[mol][1])

        old_method = False
        if old_method:
            FES[mol][2] = FES[mol][2] - np.amax(FES[mol][2])
            conts = np.arange(-vmax, 0.0, 2.0)
        else:
            # conts = np.append(np.arange(0.0,20.0,2.0), vmax)
            conts = np.arange(0.0, vmax, 2.0)
            print(conts)
        f_x = np.linspace(0.0, 4.5, 1000)  # funnel lower & upper walls
        sc = 2.5 if old else 3.0  # funnel s-cent
        print(sc)
        b = 1.5  # funnel beta-cent
        f = 0.15  # funnel wall buffer
        h = 1.7 if old else 1.2  # funnel wall width
        f_y = h * (1.0 / (1.0 + np.exp(b * (f_x - sc)))) + f

        plt.figure()
        plt.contourf(
            x,
            y,
            FES[mol][2] / 4.184,
            conts,
            cmap="CMRmap",
        )
        plt.colorbar(label="Free Energy Surface / kcal/mol")
        plt.plot(f_x, f_y, "k")
        plt.xlim(-0.2, 5.0)
        plt.ylim(-0.1, 2.0)
        plt.title("{m}  |  {}  |  FES  |  {}ns".format(folder, t, m=mol))
        plt.xlabel(x_name + " / nm")
        plt.ylabel(y_name + " / nm")
        plt.savefig(
            "monitor/{m}_{f}/{m}-{f}_FES_{}.png".format(t, f=folder, m=mol),
            bbox_inches="tight",
            dpi=300,
        )

        COL = mols.copy()
        filename = "./monitor/{m}_{f}/{m}.colvar".format(f=folder, m=mol)
        with open(filename) as f:
            lines = f.readlines()
            data = [l for l in lines if l[0] not in ("@", "#")]
            data = [[float(val) for val in line.split()] for line in data]
        COL[mol] = [[]] * len(data[0])
        for i in np.arange(len(data[0])):
            COL[mol][i] = [l[i] for l in data]
        plt.figure()
        plt.plot([x / 1000 for x in COL[mol][0]], COL[mol][1], label=mol)
        plt.axhline(0.0, color="k", linestyle="--")
        plt.axhline(4.5, color="k", linestyle="--")
        plt.legend()
        plt.xlabel("Simulation Time / ns")
        plt.ylabel(x_name + " / nm")
        plt.ylim(-0.2, 5.0)
        plt.title(
            "{m}  |  {}  |  Projection on Z - pp.proj  |  {}ns".format(folder, t, m=mol)
        )
        plt.savefig(
            "monitor/{m}_{f}/{m}-{f}_proj_{}.png".format(t, f=folder, m=mol),
            bbox_inches="tight",
            dpi=300,
        )

        plt.figure()
        plt.plot([x / 1000 for x in COL[mol][0]], COL[mol][2], label=mol)
        plt.axhline(0.0, color="k", linestyle="--")
        plt.axhline(2.0, color="k", linestyle="--")
        plt.legend()
        plt.xlabel("Simulation Time / ns")
        plt.ylabel(y_name + " / nm")
        plt.ylim(-0.1, 2.1)
        plt.title(
            "{m}  |  {}  |  Distance from Z - pp.ext  |  {}ns".format(folder, t, m=mol)
        )
        plt.savefig(
            "monitor/{m}_{f}/{m}-{f}_ext_{}.png".format(t, f=folder, m=mol),
            bbox_inches="tight",
            dpi=300,
        )

        from matplotlib.ticker import NullFormatter

        plt.gcf().clear()
        newx = np.random.randn(1000)
        newy = np.random.randn(1000)

        nullfmt = NullFormatter()  # no labels

        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left + width + 0.02

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.25]
        rect_histy = [left_h, bottom, 0.255, height]

        fig = plt.figure(figsize=(12, 8))

        axContour = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)

        # no labels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)

        # the scatter plot:
        c = axContour.contourf(
            x,
            y,
            FES[mol][2] / 4.184,
            conts,
            cmap="CMRmap",
        )
        plt.colorbar(c, ax=axHisty, label="Free Energy Surface / kcal/mol")
        axContour.plot(f_x, f_y, "k")

        fig.suptitle("{m}  |  {}  |  FES  |  {}ns".format(folder, t, m=mol))

        # now determine nice limits by hand:
        binwidth = 0.25
        xymax = max(np.max(np.abs(newx)), np.max(np.abs(newy)))
        lim = (int(xymax / binwidth) + 1) * binwidth

        axContour.set_xlim((-0.2, 5.0))
        axContour.set_ylim((-0.1, 2.0))
        axContour.set_xlabel(x_name + " / nm")
        axContour.set_ylabel(y_name + " / nm")

        bins = np.arange(-lim, lim + binwidth, binwidth)
        axHistx.hist(newx, bins=bins)
        axHisty.hist(newy, bins=bins, orientation="horizontal")

        axHistx.set_xlim(axContour.get_xlim())
        axHisty.set_ylim(axContour.get_ylim())

        fig.savefig(
            "monitor/{m}_{f}/{m}-TRIPLE.png".format(t, f=folder, m=mol),
            bbox_inches="tight",
            dpi=450,
        )
