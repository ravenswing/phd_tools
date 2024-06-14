import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import seaborn as sns
from math import ceil

from numpy.polynomial import polynomial as P

from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib import ticker

from scipy.interpolate import griddata


def heatmap(data, row_labels, col_labels, ax=None, cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=False, labeltop=False, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    # original rotation = -30
    plt.setp(ax.get_xticklabels(), rotation=0, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im,
    data=None,
    valfmt="{x:.2f}",
    textcolors=("black", "white"),
    threshold=None,
    min_cutoff=None,
    max_cutoff=None,
    **textkw,
):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.0

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if min_cutoff is not None and data[i, j] < min_cutoff:
                continue
            elif max_cutoff is not None and data[i, j] > max_cutoff:
                continue
            else:
                kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
                text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                texts.append(text)

    return texts


def rmsf_gridmap(rmsf_data, seq, grid=None, offset=None):
    """
    to plot the rmsf data - specifically for B-sheet and lamellars
        rmsf_data       output from load_data.xvg
        seq             Sequence: list of residue names/letters
        grid            dimensions of heatmap grid is displaying as such
        offset          if anti-parallel, number of cells to displace each
                        alternating row by
    """
    # adjust cmap so minimum value is perfectly white
    old_cmap = cm.get_cmap("YlGnBu", 256)
    newcolors = old_cmap(np.linspace(0, 1, 256))
    white = np.array([256 / 256, 256 / 256, 256 / 256, 1])
    newcolors[0, :] = white
    custom_cmap = ListedColormap(newcolors)
    # max_cmap_val = 12.
    max_cmap_val = 20.0

    # extract data from loaded xvg pd DataFrame
    head = rmsf_data.columns.values.tolist()
    rmsf = np.array(rmsf_data[head[1]].tolist()) * 10

    # extract molecule name from sequence
    out_name = "".join([seq[-1], seq[-2], seq[1], seq[0]])

    # plots a bar or the RMSF per residue
    if grid is None:
        print("2D plot")

    elif offset is None:
        print("Plotting for parallel arrangement")
        assert len(grid) == 2
        assert rmsf_data.shape[0] == grid[0] * grid[1]

        # reshape to grid dimensions
        rmsf = np.reshape(rmsf, (grid[0], grid[1]))

        fig, ax = plt.subplots(figsize=(11, 3.5))

        im, cbar = heatmap(
            rmsf,
            [x + 1 for x in list(np.arange(grid[0]))],
            seq,
            ax=ax,
            cmap=custom_cmap,
            cbarlabel="RMSF $\AA$",
            vmin=0.0,
            vmax=max_cmap_val,
        )

        annotate_heatmap(im, valfmt="{x:.1f}", min_cutoff=-1.0)

        for s in np.arange(len(seq)):
            im.axes.text(s, -1, seq[s], horizontalalignment="center")

        out_name += "_para"

    else:
        print("Plotting for antiparallel arrangement")
        assert len(grid) == 2
        assert isinstance(offset, int)

        seq = [" "] * offset + seq

        # establish a grid of negatives that will not be shown
        base = np.ones((grid[0], grid[1] + offset)) * -10.0
        for i in range(grid[0]):
            if i % 2 == 0:
                base[i, offset:] = rmsf[i * grid[1] : (i + 1) * grid[1]]
            else:
                base[i, :-offset] = np.flip(rmsf[i * grid[1] : (i + 1) * grid[1]])

        fig, ax = plt.subplots(figsize=(12, 5))

        im, cbar = heatmap(
            base,
            [x + 1 for x in list(np.arange(6))],
            seq,
            ax=ax,
            cmap=custom_cmap,
            cbarlabel="RMSF $\AA$",
            vmin=0.0,
            vmax=max_cmap_val,
        )

        annotate_heatmap(im, valfmt="{x:.1f}", min_cutoff=-1.0)

        for s in np.arange(len(seq)):
            im.axes.text(s, -1, seq[s], horizontalalignment="center")
            im.axes.text(s - 1, grid[0], seq[-s], horizontalalignment="center")

        out_name += "_anti"

    fig.tight_layout()

    return fig, out_name

    # fig.savefig('./RMSF_{}.png'.format(out_name), dpi=300,
    # bbox_inches='tight', transparent=True)


def bubble_plot(csv, size_scale):
    """Make bubble plot from csv"""
    ddg = pd.read_csv(csv, sep=",")
    sizes = (ddg["bonds"] + 1) * size_scale
    sns.scatterplot(
        x="weight",
        y="fs2",
        s=sizes,
        data=ddg,
    )
    sns.scatterplot(x="weight", y="fs1", s=sizes, data=ddg, legend="brief")
    plt.xlabel("Molecular Weight / Da")
    plt.ylabel("Deviation from Experimental $\Delta$G / kcal/mol")
    plt.grid(alpha=0.5, zorder=1)
    plt.savefig(csv + "_bubble.png", dpi=300, transparent=True)


def old_cv_contour(fes, pdb, axes, in_vmax, ax):
    """Plot a contour plot for 2 CVs"""
    x = fes[:, 0]
    y = fes[:, 1]
    z = fes[:, 2]

    z = np.array(z / 4.184)
    z = np.subtract(z, min(z))
    max_non_inf = np.amax(z[np.isfinite(z)])
    print("VMAX: ", max_non_inf)
    x_name, y_name = axes
    # vmax = int(ceil(max_non_inf / 2.0)) * 2 if 'REW' in name else in_vmax
    vmax = in_vmax
    # vmax = 50
    x = np.array([nm * 10 for nm in x])
    y = np.array([nm * 10 for nm in y])

    xgrid = int(np.sqrt(len(x)) - 1)
    ygrid = int(np.sqrt(len(y)) - 1)

    xi = np.linspace(min(x), max(x), xgrid)
    yi = np.linspace(min(y), max(y), ygrid)

    maxz = 0
    for ndx, v in enumerate(z):
        if np.isfinite(z[ndx]):
            if z[ndx] > maxz:
                maxz = z[ndx]

    for ndx, v in enumerate(z):
        if np.isinf(z[ndx]):
            z[ndx] = maxz

    xi, yi = np.meshgrid(xi, yi)
    print(
        x.shape,
        y.shape,
        z.shape,
        xi.shape,
        yi.shape,
    )
    zi = griddata((x, y), z, (xi, yi), method="linear")

    # iso = round(2*max_non_inf/12)/2
    conts = np.arange(0.001, vmax + 1, 2.0)
    # conts = np.arange(0.001, max(z)+1, 2.0)

    f_x = np.linspace(0.0, 45, 1000)  # funnel lower & upper walls
    sc = 30
    b = 0.15  # funnel beta-cent
    f = 1.5  # funnel wall buffer
    h = 12  # funnel wall width
    f_y = h * (1.0 / (1.0 + np.exp(b * (f_x - sc)))) + f

    #    ax = fig.add_subplot(plot_n, sharex=True, sharey=True)
    CS = ax.contourf(xi, yi, zi, conts, cmap="RdYlBu", antialiased=True)
    ax.contour(
        xi, yi, zi, conts, colors="k", linewidths=0.5, alpha=0.5, antialiased=True
    )
    ax.plot(f_x, f_y, "k")
    ax.set_xlim(-2.0, 50.0)
    ax.set_ylim(-1.0, 20.0)

    ax.grid()
    return CS


def ddg_scatter(csv, mode):
    """Make custom scatter from ddG data"""
    ddg = pd.read_csv(csv, sep=",")
    plt.figure()
    line1 = P.polyfit(x=ddg["weight"], y=ddg["fs1"], deg=1)
    f1 = P.Polynomial(line1)
    line2 = P.polyfit(x=ddg["weight"], y=ddg["fs2"], deg=1)
    f2 = P.Polynomial(line2)
    line3 = P.polyfit(
        x=pd.concat([ddg["weight"], ddg["weight"]]),
        y=pd.concat([ddg["fs1"], ddg["fs2"]]),
        deg=1,
    )
    f3 = P.Polynomial(line3)
    x = np.linspace(150, 550, 50)
    plt.hlines(y=0, xmin=150, xmax=550, colors="xkcd:green", linewidth=2.5)

    if mode == 1:
        plt.errorbar(
            x="weight",
            y="fs2",
            yerr=2.0,
            data=ddg,
            fmt="o",
            capsize=5,
            c=colours[2],
            label="FS2",
        )
        plt.errorbar(
            x="weight",
            y="fs1",
            yerr=2.0,
            data=ddg,
            fmt="o",
            capsize=5,
            c=colours[0],
            label="FS1",
        )
    if mode == 2:
        plt.scatter(
            x="weight",
            y="fs2",
            data=ddg,
            marker="D",
            s=6,
            c=colours[2],
            label="FS2",
            zorder=2,
        )
        plt.scatter(
            x="weight",
            y="fs1",
            data=ddg,
            marker="D",
            s=6,
            c=colours[0],
            label="FS1",
            zorder=3,
        )
        plt.axhspan(-2, 2, facecolor="xkcd:green", alpha=0.2, zorder=1)

    plt.plot(x, f1(x), "--", c=colours[0], alpha=0.5)
    plt.plot(x, f2(x), "--", c=colours[2], alpha=0.5)
    plt.plot(x, f3(x), "k", label="Combined Trend")
    plt.xlim([150, 550])
    plt.ylim([-1.0, 15.0])
    plt.xlabel("Molecular Weight / Da")
    plt.ylabel("Deviation from Experimental $\Delta$G / kcal/mol")
    plt.grid(alpha=0.5)
    plt.legend()
    plt.savefig(
        csv.split(".")[0] + str(mode) + "_scatter.png", dpi=300, transparent=True
    )
