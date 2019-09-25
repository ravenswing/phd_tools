"""
===============================================================================
                                GRAPHICS AND PLOTTING
===============================================================================
"""

from math import ceil
import numpy as np
from numpy.polynomial import polynomial as P
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

colours = ['#31859C',   # FS1 & BS1
           '#FFC000',   # Tunnel
           '#7030A0',   # FS2 & BS2
          ]

def hills_plot(hills_data, pdb, funnel_side, save_dir):
    """ Plot a simple line graph of HILLS file """
    plt.figure()
    plt.plot([x/1000 for x in hills_data[0]], hills_data[1], label=pdb)
    plt.legend()
    plt.title('{m}  |  {}  |  Hill Heights'.format(funnel_side, m=pdb))
    plt.savefig('{}/{m}-{f}_Heights.png'\
            .format(save_dir, f=funnel_side, m=pdb), bbox_inches='tight', dpi=300)

def diffusion_plots(colvar_data, pdb, funnel_side, num_cvs, save_dir):
    """ plots the diffusion plots for X cvs."""
    if num_cvs == 2:
        x_name = "X-axis"
        y_name = "Y-axis"
        plt.figure()
        plt.plot([x/1000 for x in colvar_data[:, 0]],
                 colvar_data[:, 1], label=pdb)
        plt.axhline(0.0, color='k', linestyle='--')
        plt.axhline(4.5, color='k', linestyle='--')
        plt.legend()
        plt.xlabel('Simulation Time / ns')
        plt.ylabel(x_name+' / nm')
        plt.ylim(-0.2, 5.0)
        plt.title('{m}  |  {}  |  Projection on Z - pp.proj '\
                .format(funnel_side, m=pdb))
        plt.savefig('{}/{m}-{f}_proj.png'\
                .format(save_dir, f=funnel_side, m=pdb), bbox_inches='tight', dpi=300)
        plt.figure()
        plt.plot([x/1000 for x in colvar_data[:, 0]],
                 colvar_data[:, 2], label=pdb)
        plt.axhline(0.0, color='k', linestyle='--')
        plt.axhline(2.0, color='k', linestyle='--')
        plt.legend()
        plt.xlabel('Simulation Time / ns')
        plt.ylabel(y_name+' / nm')
        plt.ylim(-0.1, 2.1)
        plt.title('{m}  |  {}  |  Distance from Z - pp.ext'\
                .format(funnel_side, m=pdb))
        plt.savefig('{}/{m}-{f}_ext.png'\
                .format(save_dir, f=funnel_side, m=pdb), bbox_inches='tight', dpi=300)


def two_cv_contour(fes, pdb, funnel_side, axes, in_vmax, name, save_dir, ax):
    """ Plot a contour plot for 2 CVs"""

    # NEED TO ADD THIS FUNCTIONALITY BACK IN!!!
    #with open('./monitor/vmax.dat', 'a+') as f:
    #   f.write('{} {} {}\r\n'.format(pdb, funnel_side, np.amax(fes[2])))
    fes[2] = fes[2]/4.184
    max_non_inf = np.amax(fes[2][np.isfinite(fes[2])])
    print('VMAX: ', max_non_inf)
    x_name, y_name = axes
    vmax = int(ceil(max_non_inf / 2.0)) * 2 if 'REW' in name else in_vmax
    vmax = 50

    x, y = np.meshgrid(fes[0], fes[1])

    iso = round(2*max_non_inf/12)/2

    old_method = False
    if old_method:
        fes[2] = fes[2]-np.amax(fes[2])
        conts = np.arange(-vmax, 0.0, 1.0)
    else:
        conts = np.arange(0.0001, max_non_inf, iso)

    f_x = np.linspace(0.0, 4.5, 1000) # funnel lower & upper walls
    sc = 3.0
    b = 1.5     # funnel beta-cent
    f = 0.15    # funnel wall buffer
    h = 1.2     # funnel wall width
    f_y = h*(1./(1.+np.exp(b*(f_x-sc))))+f

#    ax = fig.add_subplot(plot_n, sharex=True, sharey=True)
    CS = ax.contourf(x, y, fes[2], conts, cmap='RdYlBu', antialiased=True)
    ax.contour(x, y, fes[2], conts, colors='k', linewidths=0.5, alpha=0.5, antialiased=True)
    #ax.colorbar(CS,format='%.1f',use_gridspec=True, label='Free Energy Surface / kcal/pdb')
    if 'REW' not in name:
        ax.plot(f_x, f_y, 'k')
        #ax.set_xlim(-0.2, 5.0)
        #ax.set_ylim(-0.1, 2.0)
        #ax.set_title('{m}  |  {}'.format(funnel_side, m=pdb))
    else:
        plt.xlim(-0.2, 5.0)
        plt.ylim(-0.1, 2.0)
        #ax.title('{m}  |  {}  |  Reweighted Free Energy Surface'.format(funnel_side, m=pdb))
    ax.grid()
    #ax.set_xlabel(x_name+' / nm')
    #ax.set_ylabel(y_name+' / nm')
    #ax.savefig('{}/{m}-{f}_{}.png'\
     #       .format(save_dir, name, f=funnel_side, m=pdb), bbox_inches='tight', dpi=300)


def make_ax(fig, n):
    """ Add axis (for multiplot) """
    x = np.linspace(1, 10)
    y = x * n
    ax = fig.add_subplot(2, 3, n)
    ax.plot(x, y)
    ax.set_title('Times '+str(n))
    return fig

def bubble_plot(csv, size_scale):
    """ Make bubble plot from csv """
    ddg = pd.read_csv(csv, sep=',')
    #print(ddg.head)
    sizes = (ddg["bonds"] +1)*size_scale
    sns.scatterplot(x='weight', y='fs2', s=sizes, data=ddg,)
    sns.scatterplot(x='weight', y='fs1', s=sizes, data=ddg, legend='brief')
    plt.xlabel('Molecular Weight / Da')
    plt.ylabel('Deviation from Experimental $\Delta$G / kcal/mol')
    plt.grid(alpha=0.5, zorder=1)
    plt.savefig(csv+'_bubble.png', dpi=300, transparent=True)

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
    plt.savefig(csv.split('.')[0]+str(mode)+'_scatter.png',
                dpi=300, transparent=True)
