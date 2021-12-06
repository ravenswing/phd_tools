"""
===============================================================================
                                GRAPHICS AND PLOTTING
===============================================================================
"""

from math import ceil
import re
import numpy as np
from scipy.interpolate import griddata
from numpy.polynomial import polynomial as P
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

colours = ['#31859C',   # FS1 & BS1 = RIGHT HAND SIDE
           '#FFC000',   # Tunnel
           '#7030A0',   # FS2 & BS2 = LEFT HAND SIDE
          ]

def hills_plot(hills_data, pdb, funnel_side, save_dir):
    """ Plot a simple line graph of HILLS file """
    plt.figure()
    plt.plot([x/1000 for x in hills_data[0]], hills_data[1], label=pdb)
    plt.legend()
    plt.title('{m}  |  {}  |  Hill Heights'.format(funnel_side, m=pdb))
    plt.savefig('{}/{m}-{f}_Heights.png'\
            .format(save_dir, f=funnel_side, m=pdb), bbox_inches='tight', dpi=300)

def single_diffusion_plots(colvar_data, pdb, funnel_side, num_cvs, save_dir):
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
    x = fes[:, 0]
    y = fes[:, 1]
    z = fes[:, 2]
    

    z = np.array(z/4.184)
    z = np.subtract(z, min(z))
    max_non_inf = np.amax(z[np.isfinite(z)])
    print('VMAX: ', max_non_inf)
    x_name, y_name = axes
    #vmax = int(ceil(max_non_inf / 2.0)) * 2 if 'REW' in name else in_vmax
    vmax = in_vmax
    #vmax = 50
    x = np.array([nm*10 for nm in x])
    y = np.array([nm*10 for nm in y])
    
    xgrid = int(np.sqrt(len(x))-1)
    ygrid = int(np.sqrt(len(y))-1)
   # xgrid = 300
   # ygrid = 300
    
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
    print(x.shape, y.shape, z.shape, xi.shape, yi.shape,)
    zi = griddata((x,y), z, (xi, yi), method='linear')

    iso = round(2*max_non_inf/12)/2
    
    conts = np.arange(0.001, vmax+1, 2.0)
    #conts = np.arange(0.001, max(z)+1, 2.0)

    f_x = np.linspace(0.0, 45, 1000) # funnel lower & upper walls
    sc = 30
    b = 0.15     # funnel beta-cent
    f = 1.5    # funnel wall buffer
    h = 12     # funnel wall width
    f_y = h*(1./(1.+np.exp(b*(f_x-sc))))+f

#    ax = fig.add_subplot(plot_n, sharex=True, sharey=True)
    CS = ax.contourf(xi, yi, zi, conts, cmap='RdYlBu', antialiased=True)
    ax.contour(xi, yi, zi, conts, colors='k', linewidths=0.5, alpha=0.5,
               antialiased=True)
    if 'REW' not in name:
        ax.plot(f_x, f_y, 'k')
        ax.set_xlim(-2.0, 50.0)
        ax.set_ylim(-1.0, 20.0)
        #ax.set_title('{m}  |  {}'.format(funnel_side, m=pdb))
    else:
        print('?')
        plt.xlim(0., 50.)
        plt.ylim(0., 50.)
        ax.set_ylabel('$\mathrm{RMSD_{OUT}}$ / $\mathrm{\AA}$')
        #ax.title('{m}  |  {}  |  Reweighted Free Energy Surface'.format(funnel_side, m=pdb))
    ax.grid()
    return CS

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

def convergence(fes_dir, ts_list, ax):
    """ Plot convergence of cv """
    lin_cols = ['xkcd:light red', 'xkcd:light orange', 'xkcd:light green',
                'xkcd:light cyan', 'xkcd:ocean blue']
    init_file = '{}/fes_{}.dat'.format(fes_dir, int(ts_list[0]/10))
    conv_data = pd.concat([df[df.cv != "#!"] for df in \
                          pd.read_csv(init_file, delim_whitespace=True,
                                      names=['cv', str(ts_list[0])], skiprows=5,
                                      chunksize=1000)])
    for timestamp in ts_list[1:]:
        fes_file = '{}/fes_{}.dat'.format(fes_dir, int(timestamp/10))
        fes_data = pd.concat([df[df.cv != "#!"] for df in \
                             pd.read_csv(fes_file, delim_whitespace=True,
                                         names=['cv', str(timestamp)],
                                         skiprows=5, chunksize=1000)])
        conv_data = pd.merge(conv_data, fes_data, on='cv')
    for i in np.arange(len(ts_list)):
        ax.plot(conv_data['cv']*10, [y/4.184 for y in conv_data[str(ts_list[i])]],
                c=lin_cols[i], label=str(ts_list[i])+' ns')

    nm = re.split('/|_', fes_dir)
    rew_file = '{p}_{f}/{p}-{f}_{c}.fes'.format(p=nm[0], f=nm[1], c=nm[-1])
    rew_data = pd.read_csv(rew_file, delim_whitespace=True, names=['rx', 'ry'], skiprows=5)
    #rew_data = rew_data[rew_data.rx != "#!"]
    print(rew_file)
    print(rew_data.head())
    ax.plot(rew_data['rx']*10, [y/4.184 for y in rew_data['ry']], 'k')
    if 'proj' in fes_dir:
        ax.set_xlim([-3, 50])
        ax.set_xticks(np.arange(6)*10)
    else:
        ax.set_xlim([-1, 17])
        ax.set_xticks(np.linspace(0., 15, num=4))
    ax.set_ylim([0, 20.])
    ax.grid(alpha=0.5)


def diffusion(colv_data, pdb, cv, bounds, ax, lin_col):
    """ Plot diffusion of CV """
    colv_data['pp.'+cv] = colv_data['pp.'+cv].astype(float)
    colv_data['mean'] = colv_data['pp.'+cv].rolling(5000, center=True).mean()
    y1 = colv_data['pp.'+cv].values*10
    y2 = colv_data['mean'].values*10
    print(len(y1))
    with open('./min_proj.dat', 'a') as f:
        f.write("{}    {:5.3f}    {:5.3f}\n".format(pdb, y1[0], y2[2500]))
    x = np.arange(len(y1))*0.002
    ax.plot(x, y1, c=lin_col, alpha=0.3, lw=0.5)
    ax.plot(x, y2, c='k', alpha=1., lw=1.)
    if cv == 'proj':
        ax.set_ylim([-3., 50.])
        ax.set_yticks(np.arange(6)*10)
        ax.axhline(y=30, xmax=500., c='k', alpha=0.5, lw=1., ls='--')
        ax.axhline(y=bounds[0], xmax=500., c='k', alpha=0.5, lw=1., ls='--', zorder=20)
        ax.fill_between(x, 0, bounds[0]+bounds[1], color='k', alpha=0.1)
    else:
        ax.set_ylim([-1, 17])
        ax.set_yticks(np.linspace(0., 15, num=4))
    ax.set_xlim([0, 500.])
    ax.set_xticks(np.linspace(0., 500., num=6))
    ax.grid(alpha=0.3)

def SWISH_diffusion(pdb, _lambda, ax, data, lin_col, n, label_N, bounds, demux=False):
    N = 50 if not demux else 1
    data = data.iloc[::N, :]
    cv='proj'
 
    data['pp.'+cv] = data['pp.'+cv].astype(float) 
    y1 = data['pp.'+cv].values*10
   
    x = np.arange(len(y1))*(0.002*N) if not demux else np.arange(len(y1))
    MAX = int(np.max(x))
    ax.plot(x, y1, c=lin_col, alpha=0.3, lw=0.5)
    
    ax.set_ylim([-3, 50.])
    ax.set_yticks(np.arange(6)*10)
    ax.axhline(y=30, xmax=MAX, c='k', alpha=0.5, lw=1., ls='--')
    ax.axhline(y=bounds[0], xmax=MAX, c='k', alpha=0.5, lw=1., ls='--', zorder=20)
    ax.fill_between(x, 0, bounds[0]+bounds[1], color='k', alpha=0.1)
     
    ax.set_xlim([0, MAX])
    #ax.set_xticks(np.linspace(0., MAX, num=int(MAX/5+1)))
    ax.grid(alpha=0.3)
     
    if n % 2 == 0: ax.set_ylabel('pp.proj / $\mathrm{\AA}$', va='center', rotation='vertical', labelpad=10, fontsize=10)
    if n > label_N: ax.set_xlabel('Simulation Time / ns', ha='center', fontsize=10, labelpad=10) 


def dgdt(y, exp_value, ax, err=None):

    x = np.linspace(0., len(y)*10, len(y))

    ax.scatter(x,y, c='k', s=8, marker='D', zorder=1)
    ax.plot(x,y, c='k', zorder=2, alpha=.5)
    if err is not None:
        ax.fill_between(x, y-err, y+err, color='k', alpha=0.1)

    ax.axhline(y=exp_value, xmin=0., xmax=max(x), c='xkcd:green', ls='--')

    ax.axhspan(ymin=exp_value-2, ymax=exp_value+2, xmin=0., xmax=max(x),  facecolor='xkcd:green', alpha=0.2)
    ax.axhline(y=exp_value+2, xmin=0., xmax=max(x),  color='xkcd:green', alpha=0.2)
    ax.axhline(y=exp_value-2, xmin=0., xmax=max(x), color='xkcd:green', alpha=0.2)
    
    ax.axhspan(ymin=exp_value-3.5, ymax=exp_value+3.5, xmin=0., xmax=max(x),  facecolor='xkcd:orange', alpha=0.2)
    ax.axhline(y=exp_value+3.5, xmin=0., xmax=max(x),  color='xkcd:orange', alpha=0.2)
    ax.axhline(y=exp_value-3.5, xmin=0., xmax=max(x), color='xkcd:orange', alpha=0.2)

    ax.set_xlim([-5., max(x)+10])
    ax.set_ylim([-30., 10.0])
    ax.grid(alpha=0.3)


def ddgdt(y, exp_value, ax):
    x = np.linspace(0., len(y)*10, len(y))
    y = np.abs(np.subtract(y, exp_value))
    ax.scatter(x, y, c='k', s=8, marker='D', zorder=1)
    ax.plot(x, y, c='k', zorder=2, alpha=.5)
