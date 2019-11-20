"""
                    THIS IS ACTUALLY THE PLOTTING ONE
"""
import numpy as np
from numpy.polynomial import polynomial as P
import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import datetime
from sklearn import metrics
#import seaborn as sns

colours = ['#31859C',   # FS1 & BS1
           '#FFC000',   # Tunnel
           '#7030A0',   # FS2 & BS2
          ]

labels = ["$\mathrm{F_{RHS}}$", "$\mathrm{F_{LHS}}$"]
comet_codes = ['5alg', '5am3', '5aly', '5alp', '5aia', '5alt', '5akk',
               '5akg',] # '5alx']

def ddg_scatter(csv):
    """ Make custom scatter from ddG data """
    ddg = pd.read_csv(csv, sep=',')
    plt.figure()
    line1 = P.polyfit(x=ddg['weight'], y=ddg['fs1'], deg=1)
    f1 = P.Polynomial(line1)
    line2 = P.polyfit(x=ddg['weight'], y=ddg['fs2'], deg=1)
    f2 = P.Polynomial(line2)
    line3 = P.polyfit(x=pd.concat([ddg['weight'], ddg['weight']]),
                      y=pd.concat([ddg['fs1'], ddg['fs2']]), deg=1)
    f3 = P.Polynomial(line3)
    x = np.linspace(150, 550, 50)
    plt.hlines(y=0, xmin=150, xmax=550, colors='xkcd:green', linewidth=2.5)
    plt.scatter(x='weight', y='fs2', data=ddg,
                marker='D', s=6, c=colours[2], label=labels[1], zorder=2)
    plt.scatter(x='weight', y='fs1', data=ddg,
                marker='D', s=6, c=colours[0], label=labels[0], zorder=3)
    plt.axhspan(-2, 2, facecolor='xkcd:green', alpha=0.2, zorder=1)

    plt.plot(x, f1(x), '--', c=colours[0], alpha=0.5)
    plt.plot(x, f2(x), '--', c=colours[2], alpha=0.5)
    plt.plot(x, f3(x), 'k', label='Combined Trend')
    plt.xlim([150, 550])
    plt.ylim([-1., 15.])
    plt.xlabel('Molecular Weight / Da')
    plt.ylabel('Deviation from $\mathrm{\Delta G_{exp}}$ /'
               + ' kcal $\mathrm{mol^{-1}}$')
    plt.grid(alpha=0.5)
    plt.legend()
    plt.savefig(csv.split('.')[0]+'_scatter.png',
                dpi=300, transparent=True)


def dg_scatter(ax, df, Y):
    """ Make scatter dG """
    x = np.linspace(-20, 10, 100)
    # plot data
    ax.scatter(x='exp', y=Y, data=df,
                    marker='D', s=8, c='k', zorder=2)
    # add central line and coloured regions
    ax.plot(x, x, 'k')
    ax.fill_between(x, x+2, x-2, facecolor='xkcd:green', alpha=0.1)
    ax.plot(x, x+2, c='xkcd:green', alpha=0.2 )
    ax.plot(x, x-2, c='xkcd:green', alpha=0.2)
    ax.fill_between(x, x+2, x+3.5, facecolor='xkcd:orange', alpha=0.1)
    ax.plot(x, x+3.5, c='xkcd:orange', alpha=0.2)
    ax.fill_between(x, x-2, x-3.5, facecolor='xkcd:orange', alpha=0.1)
    ax.plot(x, x-3.5, c='xkcd:orange', alpha=0.2)
    # adjust scatter
    ax.set_xlim([-20, 5.0])
    ax.set_ylim([-20, 5.0])
    ax.grid(alpha=0.5)
    # add labels
    for i, j, n in zip(df['exp'], df[Y], df['pdb']):
        if j < i+3.7 and j > i-3.7: continue
        if j > i:
            ax.annotate(n, (i, j), xytext=(-5,0), textcoords='offset points',
                        size=8, ha='right', va="bottom")
        else:
            ax.annotate(n, (i, j), xytext=(5, 0), textcoords='offset points',
                        size=8, ha='left', va="bottom")

def split_plot(csv):
    """ Make SI 2x3 plots of dG for each binding site """
    df = pd.read_csv(csv, sep=',')
    fig = plt.figure(figsize=(10, 14))
    axes = fig.subplots(3, 2, sharex='col', sharey=True)
    fig.subplots_adjust(hspace=0.08, wspace=0.08)
    sites = ['Tunnel', 'BS1', 'BS2']
    for site in [0, 1, 2]:
        to_plot = df.loc[df['site'] == sites[site]]
        for fs in [1, 2]:
            dg_scatter(axes[site, fs-1], to_plot, 'fs'+str(fs))

            fig.text(0.31+((fs-1)*0.41), 0.9, labels[fs-1], ha='center', fontsize=18)
            fig.text(0.31+((fs-1)*0.41), 0.07,
                     '$\mathrm{\Delta G_{exp}}$ / kcal $\mathrm{mol^{-1}}$',
                     ha='center', fontsize=10)

        fig.text(0.05, 0.24+(site*0.27),
                 '$\mathrm{\Delta G_{calc}}$ / kcal $\mathrm{mol^{-1}}$',
                 va='center', rotation='vertical', fontsize=10)
    fig.savefig('FunMetaD_FSsplit_dG.png', dpi=300, transparent=True)

def quad_plot(csv, SI=False):
    """ Make manuscript and SI 2x2 plots of all methods """
    #comet_codes = ['5alg', '5am3', '5aly', '5alp', '5aia', '5alt', '5akk',
#                   '5akg',] # '5alx']
    dg_data = pd.read_csv(csv, sep=',')
    fig = plt.figure(figsize=(10, 10))
    axes = fig.subplots(2, 2, sharex='col', sharey=True)
    fig.subplots_adjust(hspace=0.12, wspace=0.08)
    df = dg_data[dg_data['pdb'].isin(comet_codes)] if SI else dg_data

    dg_scatter(axes[0, 0], df, 'fs1')
    axes[0, 0].set_title('Fun-metaD - '+labels[0], fontsize=16)
    dg_scatter(axes[0, 1], df, 'fs2')
    axes[0, 1].set_title('Fun-metaD - '+labels[1], fontsize=16)
    dg_scatter(axes[1, 0], df, 'comet')
    axes[1, 0].set_title('CoMet Path', fontsize=16)
    dg_scatter(axes[1, 1], df, 'swish')
    axes[1, 1].set_title('Fun-SWISH', fontsize=16)

    for i in [0,1]:
        # x label
        fig.text(0.31+(i*0.41), 0.07,
                 '$\mathrm{\Delta G_{exp}}$ / kcal $\mathrm{mol^{-1}}$',
                 ha='center', fontsize=10)
        # y label
        fig.text(0.05, 0.31+(i*0.41),
                 '$\mathrm{\Delta G_{calc}}$ / kcal $\mathrm{mol^{-1}}$',
                 va='center', rotation='vertical', fontsize=10)

    s = 'SI' if SI else 'manuscript'
    fig.savefig('Quad_dG_{}.png'.format(s), dpi=300, transparent=True)

def stats(y_true, y_pred):
    # R-squared
    #r2 = metrics.r2_score(y_true, y_pred)
    # RMSE
    mae = metrics.mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(metrics.mean_squared_error(y_true, y_pred))
    # Pearson r
    r = y_true.corr(y_pred, method='pearson')
    r2 = r**2
    # Kendall tau
    tau = y_true.corr(y_pred, method='kendall')

    return [r2, r, tau, mae, rmse]

def write_stats(csv):
    dg_data = pd.read_csv(csv, sep=',')
    grouped = dg_data.groupby('site')

    current = datetime.datetime.now()
    logfile = current.strftime("Statistics_%d%b-%H-%M.md")
    path = pathlib.Path('./'+logfile)
    head = ("| Site | FS |  R-squared | Pearson *r* | Kendall *tau* | MAE | RMSE | N |\n"
            "|------|----|------------|-------------|---------------|-----|------|---|\n")
    with path.open(mode='w+') as f:
        f.write("## Per Site Per FS\n")
        f.write(head)
        for site, group in grouped:
            r, c = group.shape
            exp_val = group.exp
            fs1_val = group.fs1
            stat_list1 = stats(exp_val, fs1_val)
            f.write('|{}| RHS |{}| {} |\n'
                    .format(site,
                            ' | '.join(['{:8.6f}'.format(n) for n in stat_list1]),
                            r ))
            fs2_val = group.fs2
            stat_list2 = stats(exp_val, fs2_val)
            f.write('| {} | LHS | {} | {} |\n'
                    .format(site,
                            ' | '.join(['{:8.6f}'.format(n) for n in stat_list2]),
                            r))

    head = ("| Method | R-squared | Pearson *r* | Kendall *tau* | MAE | RMSE | N |\n"
            "|--------|-----------|-------------|---------------|-----|------|---|\n")
    with path.open(mode='a') as f:
        f.write("\n## Per Methodology (all sys)\n")
        f.write(head)
        for m in ['fs1', 'fs2', 'swish', 'comet']:
            df = dg_data.dropna(subset=['exp', m])
            r, c = df.shape
            exp_val = df.exp
            m_val = df[m]
            stat_list = stats(exp_val, m_val)
            f.write('| {} | {} | {} |\n'
                    .format(m,
                            ' | '.join(['{:8.6f}'.format(n) for n in stat_list]),
                            r))

    cut = dg_data[dg_data['pdb'].isin(comet_codes)]
    with path.open(mode='a') as f:
        f.write("\n## Per Methodology (9 sys only)\n")
        f.write(head)
        for m in ['fs1', 'fs2', 'swish', 'comet']:
            df = cut.dropna(subset=['exp', m])
            r, c = df.shape
            exp_val = df.exp
            m_val = df[m]
            stat_list = stats(exp_val, m_val)
            f.write('| {} | {} | {} |\n'
                    .format(m,
                            ' | '.join(['{:8.6f}'.format(n) for n in stat_list]),
                            r))

if __name__ == "__main__":

    #ddg_scatter('ddg_data.csv')
    #split_plot('dg_values.csv')
    #quad_plot('dg_values.csv')
    #quad_plot('dg_values.csv', SI=True)

    write_stats('dg_values.csv')
