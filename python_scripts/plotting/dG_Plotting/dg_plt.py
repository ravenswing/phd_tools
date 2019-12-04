import numpy as np
from numpy.polynomial import polynomial as P
import matplotlib.pyplot as plt
import pandas as pd
#import seaborn as sns

colours = ['#31859C',   # FS1 & BS1
           '#FFC000',   # Tunnel
           '#7030A0',   # FS2 & BS2
          ]

def bubble():
    ddg = pd.read_csv('./ddg_data.csv', sep=',')
    #print(ddg.head)
    sizes = (ddg["bonds"] +1)*100
    #sns.scatterplot(x='weight', y='fs2', s=sizes, data=ddg,)
    #sns.scatterplot(x='weight', y='fs1', s=sizes, data=ddg, legend='brief')
    plt.xlabel('Molecular Weight / Da')
    plt.ylabel('Deviation from Experimental $\Delta$G / kcal/mol')
    plt.grid(alpha=0.5, zorder=1)
    plt.savefig('Test_bubbleplot.png', dpi=300, transparent=True)

def ddg_scatter(csv, mode):
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

    if mode == 1:
        plt.errorbar(x='weight', y='fs2', yerr=2.0, data=ddg,
                     fmt='o', capsize=5, c=colours[2], label='FS2')
        plt.errorbar(x='weight', y='fs1', yerr=2.0, data=ddg,
                     fmt='o', capsize=5, c=colours[0], label='FS1')
    if mode == 2:
        plt.scatter(x='weight', y='fs2', data=ddg,
                    marker='D', s=6, c=colours[2], label='FS2', zorder=2)
        plt.scatter(x='weight', y='fs1', data=ddg,
                    marker='D', s=6, c=colours[0], label='FS1', zorder=3)
        plt.axhspan(-2, 2, facecolor='xkcd:green', alpha=0.2, zorder=1)

    plt.plot(x, f1(x), '--', c=colours[0], alpha=0.5)
    plt.plot(x, f2(x), '--', c=colours[2], alpha=0.5)
    plt.plot(x, f3(x), 'k', label='Combined Trend')
    plt.xlim([150, 550])
    plt.ylim([-1., 15.])
    plt.xlabel('Molecular Weight / Da')
    plt.ylabel('Deviation from Experimental $\Delta$G / kcal/mol')
    plt.grid(alpha=0.5)
    plt.legend()
    print(csv.split('.'))
    plt.savefig(csv.split('.')[0]+str(mode)+'_scatter.png',
                dpi=300, transparent=True)

def dg_scatter(ax, df):
    """ Make custom scatter dG values """

    #swish_codes = ['5alg', '5am3', '5aly', '5alp', '5aia', '5alt', '5akk',
                   #'5akg',] # '5alx']
    #dg_data = pd.read_csv(csv, sep=',')
    #df = dg_data[dg_data['pdb'].isin(swish_codes)] if swish else dg_data 
    #fig = plt.figure(figsize=(9,4))
    #axes = fig.subplots(1, 2, sharey=True)
    #axes.ravel()
    x = np.linspace(-20,10,100)

    ax.scatter(x='exp', y='fs'+str(fs), data=df,
                    marker='D', s=6, c=colours[(fs-1)*2], label='FS'+str(fs), zorder=2)
    ax.plot(x, x, '--k')
    ax.plot(x, x+2, '--g')
    ax.plot(x, x-2, '--g')
    ax.plot(x, x+3, '--r')
    ax.plot(x, x-3, '--r')

    ax.set_ylabel('Calculated $\Delta$G / kcal mol$^{-1}$')
    ax.set_xlabel('Experimental $\Delta$G / kcal mol$^{-1}$')
    ax.set_xlim([-20, 0])
    ax.set_ylim([-20, 10])
    ax.grid(alpha=0.5)

    for i, j, n in zip(df['exp'], df['fs'+str(fs)], df['pdb']):
        if j < i+3.2 and j > i-3.2: continue
        ax.annotate(n, (i, j), xytext=(-40, 10), textcoords='offset points',
        size=8, ha='left', va="bottom", arrowprops=dict(arrowstyle='-'))


        #bbox=dict(boxstyle="round", alpha=0.1),
        #arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1))
#    fig.text(0.31+((fs-1)*0.41), 0.9, 'FS'+str(fs), ha='center', fontsize=14) 
#s = '9-8' if swish else 'ALL'
    #fig.savefig('FunnelMetaD_dG_{}.png'.format(s),
                #dpi=300, transparent=True) 

def dg_plot(csv, swish=False, select_sys=False, split=False):
    "Actually plot the damn thing"

    swish_codes = ['5alg', '5am3', '5aly', '5alp', '5aia', '5alt', '5akk',
                   '5akg',] # '5alx']
    dg_data = pd.read_csv(csv, sep=',')
    df = dg_data[dg_data['pdb'].isin(swish_codes)] if swish else dg_data
    #axes = fig.subplots(1, 2, sharey=True)



    #fig = plt.figure(figsize=(9,4))



def swish_scatter(csv, swish=True):
    """ Make custom scatter dG values """

    swish_codes = ['5alg', '5am3', '5aly', '5alp', '5aia', '5alt', '5akk',
                   '5akg',] # '5alx']
    dg_data = pd.read_csv(csv, sep=',')
    df = dg_data[dg_data['pdb'].isin(swish_codes)] if swish else dg_data
    fig = plt.figure(figsize=(4.5,4))
    ax = plt.axes()
    x = np.linspace(-20,10,100)

    ax.scatter(x='exp', y='swish', data=df,
                    marker='D', s=6, c='k', label='SWISH', zorder=2)
    ax.plot(x, x, '--k')
    ax.plot(x, x+2, '--g')
    ax.plot(x, x-2, '--g')
    ax.plot(x, x+3.5, '--r')
    ax.plot(x, x-3.5, '--r')

    ax.set_ylabel('Calculated $\Delta$G / kcal mol$^{-1}$')
    ax.set_xlabel('Experimental $\Delta$G / kcal mol$^{-1}$')
    ax.set_xlim([-20, 0])
    ax.set_ylim([-20, 10])
    ax.grid(alpha=0.5)

    for i, j, n in zip(df['exp'], df['swish'], df['pdb']):
        if i-3.5 < j < i+3.5: continue
        ax.annotate(n, (i, j), xytext=(40, 10), textcoords='offset points',
        size=8, ha='left', va="bottom", arrowprops=dict(arrowstyle='-'))
        #bbox=dict(boxstyle="round", alpha=0.1),
        #arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1)) 
    s = '9-8' if swish else 'ALL'
    fig.savefig('SWISH_dG_{}.png'.format(s),
                dpi=300, transparent=True)



def fs_split(csv, swish=True):
    """ Make custom scatter dG values """

    swish_codes = ['5alg', '5am3', '5aly', '5alp', '5aia', '5alt', '5akk',
                   '5akg',] # '5alx']
    dg_data = pd.read_csv(csv, sep=',')
    df = dg_data[dg_data['pdb'].isin(swish_codes)] if swish else dg_data
    fig = plt.figure(figsize=(10, 14))
    axes = fig.subplots(3, 2, sharey=True)
    x = np.linspace(-20, 10, 100)
    sites = ['Tunnel', 'BS1', 'BS2']
    for site in [0, 1, 2]:
        to_plot = df.loc[df['site'] == sites[site]]
        for fs in [1, 2]:
            axes[site, fs-1].scatter(x='exp', y='fs'+str(fs), data=to_plot,
                            marker='D', s=6, c=colours[(fs-1)*2], label='FS'+str(fs), zorder=2)
            axes[site, fs-1].plot(x, x, '--k')
            axes[site, fs-1].plot(x, x+2, '--g')
            axes[site, fs-1].plot(x, x-2, '--g')
            axes[site, fs-1].plot(x, x+3, '--r')
            axes[site, fs-1].plot(x, x-3, '--r')

            axes[site, fs-1].set_ylabel('Calculated $\Delta$G / kcal mol$^{-1}$')
            axes[site, fs-1].set_xlabel('Experimental $\Delta$G / kcal mol$^{-1}$')
            axes[site, fs-1].set_xlim([-20, 0])
            axes[site, fs-1].set_ylim([-20, 10])
            axes[site, fs-1].grid(alpha=0.5)

            for i, j, n in zip(to_plot['exp'], to_plot['fs'+str(fs)], to_plot['pdb']):
                #if j < i+3.2 and j > i-3.2: continue
                axes[site, fs-1].annotate(n, (i, j), xytext=(-40, 10), textcoords='offset points',
                                          size=8, ha='left', va="bottom", arrowprops=dict(arrowstyle='-'))
                #bbox=dict(boxstyle="round", alpha=0.1),
                #arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1))
            fig.text(0.31+((fs-1)*0.41), 0.9, 'FS'+str(fs), ha='center', fontsize=14)

    s = '9-8' if swish else 'ALL'
    fig.savefig('Split_dG_{}.png'.format(s),
                dpi=300, transparent=True)

#ddg_scatter('ddg_data.csv', 2)
#dg_scatter('dg_values.csv', swish=False)
#dg_scatter('dg_values.csv', swish=True)
#swish_scatter('dg_values.csv', swish=False)
#fs_split('dg_values.csv', swish=False)
#swish_scatter('dg_values.csv', swish=True)
swish_diftimes([300, 500])
