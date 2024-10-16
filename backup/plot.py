"""
===============================================================================
                      CHAR CHAR CHAR CHAR CHAR CHARTS!!!
===============================================================================

    - Analysis
    - Plotting
    - Calculations
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import plotly.graph_objects as go

sys.path.append("/home/rhys/phd_tools/SAPS")
import load

ANG = "\u212b"

# global colours
CLR = {
    "ax": "rgb(69, 69, 69)",
    "ln": "rgb(116, 62, 122)",
    "f1": "rgb(251, 218, 230)",
    "f2": "rgb(255, 241, 194)",
    "A769": ["#e76f51", "#A23216"],
    "PF739": ["#f4a261", "#994B0B"],
    "SC4": ["#e9c46a", "#A07918"],
    "MT47": ["#2a9d8f", "#1E7167"],
    "MK87": ["#264653", "#13242A"],
}


def fes(
    fes_path,
    save_path,
    units="A",
    basins=None,
    funnel=None,
    basin_lables=None,
    xlims=None,
    ylims=None,
    labels=["CV1", "CV2"],
):
    """
    funnel_parms = {'lw': 0.0,
                    'uw': 4.5,
                    'sc': 2.5,
                    'b': 1.0,
                    'f': 0.15,
                    'h': 1.5}
    """
    data, lab = load.fes(fes_path, False)
    data[2] = data[2] / 4.184
    if units == "A":
        data[0] = np.multiply(data[0], 10)
        data[1] = np.multiply(data[1], 10)
    cmax = np.amax(data[2][np.isfinite(data[2])]) + 1

    fig = go.Figure(
        data=go.Contour(
            z=data[2],
            x=data[0],  # horizontal axis
            y=data[1],  # vertical axis
            colorscale="RdYlBu",
            contours=dict(start=0, end=cmax, size=2),
            colorbar=dict(title="Free Energy (kcal/mol)", titleside="right"),
        )
    )
    # format axes
    fig.update_xaxes(
        showline=True,
        linecolor=CLR["ax"],
        title_text=labels[0],
        linewidth=0.5,
        title_standoff=20,
        ticks="outside",
        minor_ticks="outside",
    )
    fig.update_yaxes(
        showline=True,
        linecolor=CLR["ax"],
        title_text=labels[1],
        linewidth=0.5,
        title_standoff=20,
        ticks="outside",
        minor_ticks="outside",
    )
    if units == "A":
        fig.update_xaxes(dtick=5.0)
        fig.update_yaxes(dtick=2.0)
    if units == "nm":
        fig.update_xaxes(dtick=0.5)
        fig.update_yaxes(dtick=0.2)

    # format the rest of the figure
    fig.update_layout(
        height=1600,
        width=2200,
        title_text="",
        font=dict(color=CLR["ax"], family="Arial", size=32),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        showlegend=False,
    )
    if xlims:
        fig.update_layout(xaxis_range=xlims)
    if ylims:
        fig.update_layout(yaxis_range=ylims)

    if basins:
        for i, b in enumerate(basins):
            if basin_lables is None:
                fig.add_shape(
                    type="rect",
                    x0=b[0],
                    x1=b[1],
                    y0=b[2],
                    y1=b[3],
                    line_dash="dash",
                    line=dict(color=CLR["ax"], width=5),
                )
            else:
                fig.add_shape(
                    type="rect",
                    x0=b[0],
                    x1=b[1],
                    y0=b[2],
                    y1=b[3],
                    line_dash="dash",
                    label=dict(text=str(i), font=dict(size=64)),
                    line=dict(color=CLR["ax"], width=5),
                )

    fig.write_image(save_path, scale=2)

    # ax.set_xlabel(
    # ax.set_ylabel(
    # b1 = plt.Rectangle((basins[f'{lig}-bnd'][0]/10, basins[f'{lig}-bnd'][2]/10),

    # cbar.set_label('Free Energy (kcal/mol)')


def multiplot(
    hdf_path,
    DIVS,
    save_frmt,
    title_frmt,
    plot_index=0,
    labels=None,
    xlims=None,
    ylims=None,
    average=False,
):
    data = pd.read_hdf(hdf_path, key="df")

    plots = DIVS.pop(plot_index)
    print(plots)
    print(DIVS)

    for p, plot in enumerate(plots):
        # initiate plots and titles
        fig, ax = plt.subplots(
            len(DIVS[0]), len(DIVS[1]), figsize=(len(DIVS[1]) * 8, len(DIVS[0]) * 5)
        )
        plt.suptitle(title_frmt.format(plot))
        plt.subplots_adjust(top=0.95)
        # add the plots to the axes
        for i, t1 in enumerate(DIVS[0]):
            for j, t2 in enumerate(DIVS[1]):
                if plot_index == 0:
                    df = data[plot][t1][t2]
                elif plot_index == 1:
                    df = data[t1][plot][t2]
                elif plot_index == 2:
                    df = data[t1][t2][plot]
                else:
                    print("UNSUPPORTED PLOT INDEX")
                mean = pd.Series(df).rolling(1000, center=True).mean()
                stdev = pd.Series(df).rolling(1000, center=True).std()

                upr = mean.add(stdev)
                lwr = mean.sub(stdev)

                col = "xkcd:navy"

                ax[i, j].plot(
                    df.index * 0.001, mean, c=col, lw=0.8, label=f"{t1} - {t2}"
                )
                ax[i, j].fill_between(
                    df.index * 0.001, upr.values, lwr.values, alpha=0.3
                )

                if average:
                    ax[i, j].axhline(y=df.mean(), c="k", alpha=0.5, lw=1.5, ls="--")

                # ax[i, j].set_title(f'{t1} - {t2}')
                ax[i, j].legend()
                if xlims:
                    ax[i, j].set_xlim(xlims)
                if ylims:
                    ax[i, j].set_ylim(ylims)
                if labels and i == 1:
                    ax[i, j].set_xlabel(labels[0])
                if labels and j == 0:
                    ax[i, j].set_ylabel(labels[1])
        fig.savefig(save_frmt.format(plot), dpi=450, bbox_inches="tight")
        plt.close()


def all_columns(
    hdf_path,
    save_frmt,
    title="Highlight:",
    labels=None,
    xlims=None,
    ylims=None,
    hline=False,
):
    df = pd.read_hdf(hdf_path, key="df")
    df.columns = df.columns.map("_".join)
    for col, data in df.items():
        fig, ax = plt.subplots(figsize=(12, 6.75))
        mean = data.rolling(1000, center=True).mean()
        stdev = data.rolling(1000, center=True).std()

        upr = mean.add(stdev)
        lwr = mean.sub(stdev)

        CLR = "xkcd:navy"

        ax.plot(df.index * 0.001, mean, c=CLR, lw=0.8)
        ax.fill_between(df.index * 0.001, upr.values, lwr.values, alpha=0.2)

        if hline:
            if hline == "average":
                ax.axhline(y=data.mean(), c="k", alpha=0.5, lw=1.5, ls="--")
            else:
                ax.axhline(y=hline, c="k", alpha=0.5, lw=1.5, ls="--")
        if xlims:
            ax.set_xlim(xlims)
        if ylims:
            ax.set_ylim(ylims)
        if labels:
            ax.set_xlabel(labels[0])
            ax.set_ylabel(labels[1])
        ax.set_title(f"{title} {' '.join(col.split('_'))}")
        fig.savefig(save_frmt.format(col), dpi=450, bbox_inches="tight")
        plt.close()


def cvs(
    colvar_path, save_path, cvs, units="A", title="CV Diffusion", xlims=None, ylims=None
):
    colvar = load.colvar(f"{colvar_path}")
    N = len(list(cvs.keys()))
    fig, ax = plt.subplots(1, N, figsize=(8 * N + 2, 6), layout="constrained")
    for i, cv in enumerate(cvs):
        ax[i].scatter(
            colvar.time.multiply(0.001),
            colvar[cv].multiply(10),
            c="#089682",
            s=8,
            alpha=0.4,
        )
        ax[i].set_xlabel("Time (ns)")
        ax[i].set_ylabel(f"{cvs[cv]}")
    fig.suptitle(f"{title}", fontsize=16)

    fig.savefig(f"{save_path}", bbox_inches="tight", dpi=300)
    plt.close()


def diffusion(DIVS, path_frmt, save_frmt, shape, cvs):
    for cv in cvs:
        fig, ax = plt.subplots(
            shape[0],
            shape[1],
            figsize=(shape[1] * 8, shape[0] * 5),
            sharey=True,
            sharex=True,
        )
        ax = ax.ravel()
        fig.tight_layout(h_pad=4)
        plt.suptitle(f"CV Diffusion - {cvs[cv][0]}")
        plt.subplots_adjust(top=0.94)
        for i, p in enumerate(DIVS):
            colvar = load.colvar(path_frmt.format(p=p))
            ax[i].scatter(colvar.time.multiply(0.001), colvar[cv])
            ax[i].set_title(" ".join(p))
            ax[i].set_xlabel("Simulation Time / ns")
            ax[i].set_ylabel(" / ".join(cvs[cv]))
        fig.savefig(save_frmt.format(cvs[cv][0]), dpi=300, bbox_inches="tight")
        plt.close()


def hills(
    hills_path, save_path, units="A", title="CV Diffusion", xlims=None, ylims=None
):
    hills = load.hills(f"{hills_path}")
    fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")
    ax.plot([x / 1000 for x in hills[0]], hills[1], c="#089682")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Hill Height")
    ax.set_title(f"{title}", fontsize=16)
    fig.savefig(f"{save_path}", bbox_inches="tight", dpi=300)
    plt.close()


def hills_multi(DIVS, path_frmt, save_name, shape):
    fig, ax = plt.subplots(
        shape[0],
        shape[1],
        figsize=(shape[1] * 8, shape[0] * 5),
        sharey=True,
        sharex=True,
    )
    ax = ax.ravel()
    fig.tight_layout(h_pad=4)
    plt.suptitle("Hill Heights")
    plt.subplots_adjust(top=0.94)
    for i, p in enumerate(DIVS):
        data = load.hills(path_frmt.format(p=p))
        ax[i].plot([x / 1000 for x in data[0]], data[1], c="#089682")
        ax[i].set_title(" ".join(p))
        ax[i].set_xlabel("Simulation Time / ns")
        ax[i].set_ylabel("Hills Heights")
    fig.savefig(save_name, dpi=300, bbox_inches="tight")
    plt.close()


def xvg_line(xvg_data, ax, col, line="solid", rollavg=True, label=None):
    head = xvg_data.columns.values.tolist()
    xvg_data["mean"] = xvg_data[head[1]].rolling(500, center=True).mean()
    y1 = xvg_data[head[1]].values * 10
    y2 = xvg_data["mean"].values * 10
    x = xvg_data[head[0]] / 1000
    ax.plot(x, y1, c=col, alpha=0.3, lw=0.5)
    ax.plot(x, y2, c=col, alpha=1.0, lw=1.0, ls=line, label=label)
    ax.grid(alpha=0.3)


def convergence(fes_dir, ts_list, ax):
    """Plot convergence of cv"""
    lin_cols = [
        "xkcd:light red",
        "xkcd:light orange",
        "xkcd:light green",
        "xkcd:light cyan",
        "xkcd:ocean blue",
    ]
    init_file = "{}/fes_{}.dat".format(fes_dir, int(ts_list[0] / 10))
    conv_data = pd.concat(
        [
            df[df.cv != "#!"]
            for df in pd.read_csv(
                init_file,
                delim_whitespace=True,
                names=["cv", str(ts_list[0])],
                skiprows=5,
                chunksize=1000,
            )
        ]
    )
    for timestamp in ts_list[1:]:
        fes_file = "{}/fes_{}.dat".format(fes_dir, int(timestamp / 10))
        fes_data = pd.concat(
            [
                df[df.cv != "#!"]
                for df in pd.read_csv(
                    fes_file,
                    delim_whitespace=True,
                    names=["cv", str(timestamp)],
                    skiprows=5,
                    chunksize=1000,
                )
            ]
        )
        conv_data = pd.merge(conv_data, fes_data, on="cv")
    for i in np.arange(len(ts_list)):
        ax.plot(
            conv_data["cv"],
            [y / 4.184 for y in conv_data[str(ts_list[i])]],
            c=lin_cols[i],
            label=str(ts_list[i]) + " ns",
        )

    nm = re.split("/|_", fes_dir)
    rew_file = "{p}_{f}/{p}-{f}_{c}.fes".format(p=nm[0], f=nm[1], c=nm[-1])
    rew_data = pd.read_csv(
        rew_file, delim_whitespace=True, names=["rx", "ry"], skiprows=5
    )
    ax.plot(rew_data["rx"], [y / 4.184 for y in rew_data["ry"]], "k")
    if "proj" in fes_dir:
        ax.set_xlim([-0.3, 5.0])
        ax.set_xticks(np.arange(6))
    else:
        ax.set_xlim([-0.1, 1.7])
        ax.set_xticks(np.linspace(0.0, 1.5, num=4))
    ax.set_ylim([0, 20.0])
    ax.grid(alpha=0.5)


def dgdt(y, exp_value, ax):
    x = np.linspace(0.0, len(y) * 10, len(y))

    ax.scatter(x, y, c="k", s=8, marker="D", zorder=1)
    ax.plot(x, y, c="k", zorder=2, alpha=0.5)

    ax.axhline(y=exp_value, xmin=0.0, xmax=max(x), c="xkcd:green", ls="--")

    ax.axhspan(
        ymin=exp_value - 2,
        ymax=exp_value + 2,
        xmin=0.0,
        xmax=max(x),
        facecolor="xkcd:green",
        alpha=0.2,
    )
    ax.axhline(y=exp_value + 2, xmin=0.0, xmax=max(x), color="xkcd:green", alpha=0.2)
    ax.axhline(y=exp_value - 2, xmin=0.0, xmax=max(x), color="xkcd:green", alpha=0.2)
    ax.axhspan(
        ymin=exp_value - 3.5,
        ymax=exp_value + 3.5,
        xmin=0.0,
        xmax=max(x),
        facecolor="xkcd:orange",
        alpha=0.2,
    )
    ax.axhline(y=exp_value + 3.5, xmin=0.0, xmax=max(x), color="xkcd:orange", alpha=0.2)
    ax.axhline(y=exp_value - 3.5, xmin=0.0, xmax=max(x), color="xkcd:orange", alpha=0.2)

    ax.set_xlim([-5.0, max(x) + 10])
    ax.set_ylim([-25.0, 5.0])
    ax.grid(alpha=0.3)
