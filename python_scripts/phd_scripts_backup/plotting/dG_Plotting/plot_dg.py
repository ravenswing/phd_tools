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
from scipy.stats import linregress
from sklearn.linear_model import LinearRegression
import seaborn as sns

colours = ['#31859C',   # FS1 & BS1
           '#FFC000',   # Tunnel
           '#7030A0',   # FS2 & BS2
           ]

labels = ["$\mathrm{F_{RHS}}$", "$\mathrm{F_{LHS}}$"]
comet_codes = ['5alg', '5am3', '5aly', '5alp', '5aia', '5alt', '5akk',
               '5akg', '5alx']


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
    ax.plot(x, x+2, c='xkcd:green', alpha=0.2)
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
            ax.annotate(n, (i, j), xytext=(-5, 0), textcoords='offset points',
                        size=8, ha='right', va="bottom")
        else:
            ax.annotate(n, (i, j), xytext=(5, 0), textcoords='offset points',
                        size=8, ha='left', va="bottom")


def split_plot(csv):
    """ Make SI 2x3 plots of dG for each binding site """
    df = pd.read_csv(csv, sep=',')
    fig = plt.figure(figsize=(6, 8.4))
    axes = fig.subplots(3, 2, sharex='col', sharey=True)
    fig.subplots_adjust(hspace=0.1, wspace=0.1)
    sites = ['Tunnel', 'BS1', 'BS2']
    for site in [0, 1, 2]:
        to_plot = df.loc[df['site'] == sites[site]]
        for fs in [1, 2]:
            dg_scatter(axes[site, fs-1], to_plot, 'fs'+str(fs))

            fig.text(0.31+((fs-1)*0.41), 0.9, labels[fs-1], ha='center',
                     fontsize=18)
            fig.text(0.31+((fs-1)*0.41), 0.05,
                     '$\mathrm{\Delta G_{exp}}$ / kcal $\mathrm{mol^{-1}}$',
                     ha='center', fontsize=10)

        fig.text(0.02, 0.24+(site*0.27),
                 '$\mathrm{\Delta G_{calc}}$ / kcal $\mathrm{mol^{-1}}$',
                 va='center', rotation='vertical', fontsize=10)
    fig.savefig('FunMetaD_FSsplit_dG.png', dpi=300, transparent=True)


def quad_plot(csv, SI=False):
    """ Make manuscript and SI 2x2 plots of all methods """
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
    axes[1, 0].set_title('COMet Path', fontsize=16)
    dg_scatter(axes[1, 1], df, 'swish')
    axes[1, 1].set_title('Fun-SWISH', fontsize=16)

    for i in [0, 1]:
        # x label
        fig.text(0.31+(i*0.41), 0.07,
                 '$\mathrm{\Delta G_{exp}}$ / kcal $\mathrm{mol^{-1}}$',
                 ha='center', fontsize=10)
        # y label
        fig.text(0.05, 0.31+(i*0.41),
                 '$\mathrm{\Delta G_{calc}}$ / kcal $\mathrm{mol^{-1}}$',
                 va='center', rotation='vertical', fontsize=10)

    s = 'SI' if SI else 'manuscript'
    fig.savefig('./Quad_dG_{}.png'.format(s),
                dpi=300, transparent=True)


def ind_scatterplot(csv, method):
    """ Make indvidual scatter plots of all methods """
    df = pd.read_csv(csv, sep=',')
    fig = plt.figure(figsize=(3, 3))
    axes = fig.subplots(1)

    dg_scatter(axes, df, method) 
    axes.set_xlabel('$\mathrm{\Delta G_{exp}}$ / kcal $\mathrm{mol^{-1}}$',
                    fontsize=11)

    axes.set_ylabel('$\mathrm{\Delta G_{calc}}$ / kcal $\mathrm{mol^{-1}}$',
                    fontsize=11)

    fig.savefig('./Ind_dG_{}.png'.format(method),
                dpi=300, transparent=True, bbox_inches='tight')


def swish_diftimes(timestamps, exclusions=None):
    """ Make custom scatter dG values """
    fig = plt.figure(figsize=(4.5, 4))
    ax = plt.axes()
    x = np.linspace(-20, 10, 100)
    # define marker style and size for each series
    markers = {300: ['D', 8], 500: ['x', 28], 1000: ['+', 38]}
    # load in base data set
    df = pd.read_csv('dg_values.csv', sep=',')
    for t in timestamps:
        # load new data
        dg_data = df.copy(deep=True) if t == 300 else \
                pd.read_csv('SWISH_dG_{}ns.csv'.format(t), sep=',')
        # replace data values
        if t > 301:
            for p in dg_data.pdb.values:
                #print(p)
                #print(df.loc[df['pdb'] == p, ['swish']])
                #print(dg_data)
                df.loc[df['pdb'] == p, ['swish']] = dg_data.loc[dg_data['pdb'] == p, ['r1']].values[0][0]
        # filter out exclusions 
        df_ex = df[~df['pdb'].isin(exclusions[t])]\
            if exclusions is not None else df
        # calculate lin. regression R-squared
        r2 = linregress(df_ex['exp'], df_ex['swish'])[2]**2
        # plot
        ax.scatter(x='exp', y='swish', data=df_ex,
                   marker=markers[t][0],
                   s=markers[t][1],
                   c='k',
                   lw=1,
                   label='{}: $\mathrm{{R^2}}$ = {:3.2f}'.format(t, r2),
                   zorder=2)
        ax.legend()

    #ax.scatter(x='exp', y='swish', data=df1_ex,
               #marker=markers[t],
               #s=8 if t == 300 else 28,
               #c='k',
               #label='{}: $\mathrm{{R^2}}$ = {:3.2f}'.format(t, r2),
               #zorder=2)
    ax.plot(x, x, 'k')
    ax.fill_between(x, x+2, x-2, facecolor='xkcd:green', alpha=0.1)
    ax.plot(x, x+2, c='xkcd:green', alpha=0.2)
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
    ax.set_ylabel('Calculated $\Delta$G / kcal mol$^{-1}$')
    ax.set_xlabel('Experimental $\Delta$G / kcal mol$^{-1}$')
    s = '_ALL' if exclusions is None else '_EXC'
    fig.savefig('SWISH_Times_dG{}.png'.format(s),
                dpi=300, transparent=True, bbox_inches='tight')


def stats(y_true, y_pred):
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

    combnd = pd.concat([grouped.get_group('BS1'), grouped.get_group('BS2')])
    print(combnd)
    X = combnd['exp'].to_numpy().reshape(-1, 1)
    LS = LinearRegression()
    LS.fit(X, combnd.fs1)
    print('COMBINED: ', LS.score(X, combnd.fs1))

    current = datetime.datetime.now()
    logfile = current.strftime("Statistics_%d%b-%H-%M.md")
    path = pathlib.Path('./'+logfile)
    head = ("| Site | FS |  R-squared | Pearson *r* | Kendall *tau* | MAE | RMSE | N |\n"
            "|------|----|------------|-------------|---------------|-----|------|---|\n")
    with path.open(mode='w+') as f:
        f.write("## Per Site Per FS\n")
        f.write(head)
        for site, group in grouped:
            LS = LinearRegression()
            X = group['exp'].to_numpy().reshape(-1, 1)
            LS.fit(X, group.fs1)
            print(group.site.to_numpy()[0], LS.score(X, group.fs1))
            r, c = group.shape
            exp_val = group.exp
            fs1_val = group.fs1
            stat_list1 = stats(exp_val, fs1_val)
            f.write('|{}| RHS |{}| {} |\n'
                    .format(site,
                            ' | '.join(['{:8.6f}'.format(n)
                                       for n in stat_list1]), r))
            fs2_val = group.fs2
            stat_list2 = stats(exp_val, fs2_val)
            f.write('| {} | LHS | {} | {} |\n'
                    .format(site,
                            ' | '.join(['{:8.6f}'.format(n)
                                       for n in stat_list2]), r))

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
                            ' | '.join(['{:8.6f}'.format(n)
                                       for n in stat_list]), r))

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
                            ' | '.join(['{:8.6f}'.format(n)
                                       for n in stat_list]), r))


def logP(csv):
    df = pd.read_csv(csv, sep=',')

    x = [0.05, 0.00, -0.05, -0.10, -0.15, -0.20]
    A = df.filter(regex=("r.*")).apply(lambda y: linregress(x, y), axis=1)
    df['m'] = [s[0] for s in A]
    df['err'] = [s[4] for s in A]
    df['R-sq'] = [s[2]**2 for s in A]
    df['w'] = df['err'].apply(lambda e: 1/e)
    print(df)

    X = df['logP'].to_numpy().reshape(-1, 1)

    WLS = LinearRegression()
    WLS.fit(X, df['m'], sample_weight=df['R-sq'])
    R = WLS.score(X, df['m'], sample_weight=df['R-sq'])
    print(WLS.intercept_, WLS.coef_)

    fig = plt.figure(figsize=(4.5, 4))
    ax = plt.axes()
    x = np.linspace(0.5, 5., 100)
    ax.grid(alpha=0.5)
    ax.errorbar(df['logP'], df['m'], yerr=df['err'], fmt='D', ms=6,
                c='xkcd:dark cyan')
    ax.plot(x, WLS.coef_*x + WLS.intercept_, c='xkcd:navy')
    ax.axhline(y=0., xmin=0., xmax=1, c='k', ls='--', lw=1)
    ax.set_ylabel('$\mathrm{H_{diff}}$')
    ax.set_xlabel('log P')
    fig.text(0.3, 0.8, '$\mathrm{{ R^2 }}$ = {:3.2f}'.format(R),
             ha='center', fontsize=10)
    fig.savefig('logP.png',
                dpi=300, transparent=True, bbox_inches='tight')
    plt.close(fig)

    f, ax = plt.subplots(figsize=(4.5, 4))

    # Draw the two density plots
    low_clust = df[df['logP'] < 2.5]
    high_clust = df[df['logP'] >= 2.5]

    ax.axhline(y=0., xmin=0., xmax=1, c='k', ls='--', lw=1, alpha=.8)
    ax = sns.kdeplot(low_clust['logP'], low_clust['m'],
                     cmap="Reds", shade=True, shade_lowest=False)
    ax = sns.kdeplot(high_clust['logP'], high_clust['m'],
                     cmap="Blues", shade=True, shade_lowest=False)
    ax.set_ylabel('$\mathrm{H_{diff}}$')
    f.savefig('logP_v2.png',
              dpi=300, transparent=True, bbox_inches='tight')


if __name__ == "__main__":
    #ddg_scatter('ddg_data.csv')
    split_plot('dg_values.csv')
    #quad_plot('dg_values.csv')
    #quad_plot('dg_valuesFINAL.csv')
    #quad_plot('dg_valuesFINAL.csv', SI=True)
    
    methods = ['fs1', 'fs2', 'comet', 'swish']
    for m in methods:
        ind_scatterplot('dg_valuesFINAL.csv', m)
    
    #swish_diftimes([300, 500, 1000])
    #swish_diftimes([300, 500, 1000],
                   #{300: ['5alo', '5aly', '5alh', '5alg', '5am0', '5am3'],
                    #500: ['5aly', '5alh', '5am0', '5am3'],
                    #1000: ['5alh'],
                    #})
    #logP('logP_data.csv')

    #write_stats('dg_valuesFINAL.csv')
    #write_stats('dg_values500.csv')
