"""
===============================================================================
                                    PLOTTING RUN

===============================================================================
"""

import glob
import subprocess
import sys
import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from math import floor

# local files
import graphics as gr
import load_data as load

colours = ['#31859C',   # FS1 & BS1 = RIGHT HAND SIDE 
           '#FFC000',   # Tunnel
           '#7030A0',   # FS2 & BS2 = LEFT HAND SIDE
          ]
fun_labels = {'FS1':"$\mathrm{F_{RHS}}$",
              'Tunnel':"$\mathrm{Tunnel}$",
              'FS2':"$\mathrm{F_{LHS}}$",
             }

c2p = {'5am4':0, '5aly':0, '5akh':0, '5am0':0, '5alx':0, '5aia':0, 
       '5am1':1, '5am3':1, '5alg':1, '5alp':1, '5alh':1, '5ai5':1,  
       '5alo':2, '5alt':2, '5akg':2, '5akk':2, '5ai0':2, '5ak6':2,
      }

BOUND_STATES = {'tunnel': [[10.46, 1.56], [8.96, 0.54]],
                'RHS': [[15.45, 0.94], [9.26, 1.44]],
                'LHS': [[6.86, 3.81], [10.76, 2.38]]
                }

PLOT_ALL = False
PDB = "5akk"
FS = 1
#FS = None

if PLOT_ALL:
    data_dirs = glob.glob("5a*_FS*")
elif FS is None:
    data_dirs = [PDB+"_FS1",PDB+"_FS2"]
elif FS <= 2:
    data_dirs = ["{}_FS{}".format(PDB,FS)]
else:
    print("ERROR: not sure what to plot")
    sys.exit()

def all_plots(dir_list):
    ''' plot ALL THE THINGS '''
    for directory in dir_list:
        # extract system info.
        pdb, fs = directory.split("_")
        # make directory for plots
        subprocess.call("mkdir Figures/{}/".format(directory), shell=True)
        # Hills plots
        hills_data = load.hills("{}/{}-{}.hills".format(directory, pdb, fs))
        gr.hills_plot(hills_data, pdb, fs, "Figures/{}".format(directory))
        # Diffusion plots
        diff_data = load.colvar("{}/{}-{}_OLD.colvar".format(directory, pdb, fs), 'as_numpy')
        gr.diffusion_plots(diff_data, pdb, fs, 2, "Figures/{}".format(directory))
        # OLD FES plots
        fes_data, axis_labels = load.fes("{}/{}-{}_OLD.fes".format(directory, pdb, fs), False)
        gr.two_cv_contour(fes_data, pdb, fs, axis_labels, 30, 'OLD_FES', "Figures/{}".format(directory))
        # REW FES plots
        fes_data, axis_labels = load.fes("{}/{}-{}_REW.fes".format(directory, pdb, fs), True)
        gr.two_cv_contour(fes_data, pdb, fs, ['RMSD-IN','RMSD-OUT'], 30, 'REW_FES', "Figures/{}".format(directory))

def fes_multiplot(pdb_list, cbar_max):
    fig = plt.figure(figsize=(14,14))
    axes = fig.subplots(3,2,sharex=True,sharey=True)
    n = 0
    for pdb in pdb_list:
        f = 0
        for fs in ['FS1','FS2']:
            directory = '{}_{}'.format(pdb, fs)
            plot_num = int(''.join(['2',str(len(pdb_list)),str(n)]))
            print(plot_num)
            fes_data, axis_labels = load.fes("{}/{}-{}_OLD.fes".format(directory, pdb, fs), False)   
            cmap = gr.two_cv_contour(fes_data, pdb, fs, axis_labels, cbar_max, 'OLD_FES', "Figures/{}".format(directory), axes[n,f])
            f+=1
        n += 1

    plt.subplots_adjust(left=0.15)
    # Y LABELS
    A = np.arange(len(pdb_list))
    for i in A:
        fig.text(0.1, 0.24+(A[-i-1]*0.26), axis_labels[1]+" / $\mathrm{\AA}$", va='center', rotation='vertical', fontsize=10)
        fig.text(0.06, 0.23+(A[-i-1]*0.26), pdb_list[i], ha='center', fontsize=14)
    # X LABELS
    for i in [0,1]:
        fig.text(0.34+(i*0.39), 0.065, axis_labels[0]+" / $\mathrm{\AA}$", ha='center', fontsize=10)
        fig.text(0.34+(i*0.39), 0.9, fun_labels['FS'+str(i+1)], ha='center', fontsize=18, color=colours[i*2])
    fig.subplots_adjust(hspace=0.05, wspace=0.05)
    fig.subplots_adjust(right=0.915)
    cax = plt.axes([0.93, 0.11, 0.01, 0.77])
    cbar = plt.colorbar(cmap,cax=cax, aspect=10, ticks=np.arange(0.,cbar_max+1,2.0))
    cbar.set_label("Free Energy / kcal $\mathrm{mol^{-1}}$", fontsize=10) 
    fig.savefig('Figures/FES_multi_'+'_'.join(pdb_list)+'.png', dpi=300, bbox_inches='tight')

def fes_highlight(pdb):
    """ Plot both Raw and Reweighted data in once muliplot for one system"""
    fig = plt.figure(figsize=(10, 8))
    axes = fig.subplots(2, 2, sharex='col', sharey='col')
    cbar_max = 20
    fes_list = ['OLD', 'REW']

    n = 0
    for fes_type in fes_list:
        f = 0
        is_rew = True if fes_type == 'REW' else False
        for fs in ['FS1', 'FS2']:
            directory = '{}_{}'.format(pdb, fs)
            plot_num = int(''.join(['2', str(len(fes_list)), str(n)]))
            print(plot_num)
            fes_data = load.fes_simple("{}/{}-{}_{}.fes".format(directory, pdb, fs, fes_type), is_rew)
            axis_labels = ['pp.proj', 'pp.ext'] if not is_rew else ['RMSD OUT', 'RMSD IN']
            cmap = gr.two_cv_contour(fes_data, pdb, fs, axis_labels, cbar_max, 'OLD_{}'.format(fes_type), "Figures/{}".format(directory), axes[f,n])
            f += 1
        n += 1 
    fig.subplots_adjust(hspace=0.05, wspace=0.17, left=0.15, right=0.915)
    labels = {'OLD':['Raw Data', 'pp.proj / $\mathrm{\AA}$'],
              'REW':['Reweighted', '$\mathrm{RMSD_{IN}}$ / $\mathrm{\AA}$'],
             }
    # Y LABELS
    for i in [0, 1]:
        fig.text(0.04, 0.3+(i*0.38), 'FS'+str(i+1), ha='center', fontsize=14)
        fig.text(0.08, 0.31+(i*0.38), 'pp.ext / $\mathrm{\AA}$', va='center', rotation='vertical', fontsize=10)

    # X LABELS
    for i in np.arange(len(fes_list)):
        fig.text(0.33+(i*0.41), 0.9, labels[fes_list[i]][0], ha='center', fontsize=14)
        fig.text(0.33+(i*0.41), 0.065, labels[fes_list[i]][1], ha='center', fontsize=10)

    cax = plt.axes([0.93, 0.11, 0.01, 0.77])
    cbar = plt.colorbar(cmap, cax=cax, aspect=10, ticks=np.arange(0., cbar_max+1, 2.0))
    cbar.set_label('Free Energy / kcal $\mathrm{mol^{-1}}$', fontsize=10)
    plt.show()
    fig.savefig('Figures/FES_highlight_'+pdb+'.png', dpi=300, bbox_inches='tight')

def conv_multiplot(pdb_list, mode, cv='proj'):
    """ Plot 4x2 convergence plots (one CV) for paper"""
    if mode == 'all':
        fig = plt.figure(figsize=(8, 21))
        axes = fig.subplots(4, 2, sharex=True, sharey=True)
        n = 0
        for pdb in pdb_list:
            f = 0
            for fs in ['FS1', 'FS2']:
                lines = [300, 350, 400, 450, 500] if 'kk' not in pdb else [300, 350, 400, 450, 490]
                gr.convergence('{}_{}/conv_{}'.format(pdb, fs, cv), lines, axes[n, f])
                f += 1
            n += 1
        plt.subplots_adjust(left=0.1)
        A = np.arange(len(pdb_list))
        for i in A:
            fig.text(0.04, 0.2+(A[-i-1]*0.2), pdb_list[i], ha='center', fontsize=14)
        fig.savefig('Figures/CONV_{}_multi_{}.png'.format(cv, '_'.join(pdb_list)), dpi=300, bbox_inches='tight')
        #fig.savefig('Figures/convergence_plot.png', dpi=300, bbox_inches='tight')

    elif mode == 'cut':
        cv_list = ['proj', 'ext']
        fig = plt.figure(figsize=(10, 12.5))
        axes = fig.subplots(6, 4, sharex='col', sharey=True)
        n = 0
        for pdb in pdb_list:
            f = 0
            for cv in cv_list:
                for fs in ['FS1', 'FS2']:
                    lines = [300, 350, 400, 450, 500] if 'kk' not in pdb else [300, 350, 400, 450, 490]
                    gr.convergence('{}_{}/conv_{}'.format(pdb, fs, cv), lines, axes[n, f])
                    f += 1
            n += 1

        fig.subplots_adjust(left=0.15, hspace=0.2, wspace=0.1)
        A = np.arange(len(pdb_list))
        # Y LABELS
        for i in A:
            fig.text(0.06, 0.16+(A[-i-1]*0.132), pdb_list[i], ha='center', fontsize=14, color=colours[c2p[pdb_list[i]]])
            fig.text(0.1, 0.16+(A[i]*0.132), '$\Delta$G / kcal $\mathrm{mol^{-1}}$', va='center', rotation='vertical', fontsize=10)
        # X LABELS
        for i in np.arange(len(cv_list)):
            fig.text(0.33+(i*0.38), 0.92, 'pp.'+cv_list[i], ha='center', fontsize=14)
            for j in [0, 1]:
                n = 2*i+j
                fig.text(0.24+(n*0.19), 0.065, 'pp.'+cv_list[i]+' / $\mathrm{\AA}$', ha='center', fontsize=10)
                fig.text(0.24+(n*0.19), 0.89, fun_labels['FS'+str(j+1)], ha='center', fontsize=18, color=colours[j*2])

        fig.legend([str(x)+' ns' for x in lines] + ['Reweight'], loc='lower center', ncol=6, frameon=False)
        fig.savefig('Figures/CONV_paper_'+'_'.join(pdb_list)+'.png', dpi=300, bbox_inches='tight')
    else:
        print('Mode Undefined')

def dif_multiplot(pdb_list):
    """ Plot 4x3 diffusion plots (not used)"""
    cv_list = ['proj', 'ext']
    fig = plt.figure(figsize=(10, 12.5))
    axes = fig.subplots(3, 4, sharex=True, sharey='col')
    n = 0
    for pdb in pdb_list:
        f = 0
        for cv in cv_list:
            for fs in ['FS1', 'FS2']:
                data = load.colvar('{p}_{f}/{p}-{f}_OLD.colvar'.format(p=pdb, f=fs), 'as_pandas')
                gr.diffusion(data, cv, axes[n, f])
                f += 1
        n += 1
    c2p = {'5aly':0, '5akk':2, '5alp':1, '5ai0':2, '5ai5':1, '5alt':2}
    fig.subplots_adjust(left=0.15, hspace=0.2, wspace=0.3)
    A = np.arange(len(pdb_list))
    for i in A:
        fig.text(0.06, 0.225+(A[-i-1]*0.275), pdb_list[i], ha='center', fontsize=14, color=colours[c2p[pdb_list[i]]])
        fig.text(0.1, 0.225+(A[i]*0.275), '$\Delta$G / kcal $\mathrm{mol^{-1}}$', va='center', rotation='vertical', fontsize=10)
    for i in np.arange(len(cv_list)):
        fig.text(0.33+(i*0.38), 0.92, 'pp.'+cv_list[i], ha='center', fontsize=14)
        for j in [0, 1]:
            n = 2*i+j
            fig.text(0.24+(n*0.19), 0.065, 'Simulation Time / ns', ha='center', fontsize=10)
            fig.text(0.24+(n*0.19), 0.89, 'FS'+str(j+1), ha='center', fontsize=14)
    fig.savefig('Figures/DIFF_paper.png', dpi=300, bbox_inches='tight')


def bound_check(x, bound):
    """ Establish bound, unbound of in middle """
    # check bounds are correct
    assert len(bound) == 2
    # calculate upper bound limit (avg. + std. dev.)
    threshold = bound[0] + bound[1]
    # Value of 1 = bound
    if x < threshold:
        return 1
    # Value of 2 = un-bound
    elif x > 30:
        return 2
    # Value of 0 = in the middle
    else:
        return 0


def identify_recross(input_col, bounds, max_t=500):
    """ Count the number of recrossings """

    # load in data and convert types (from str)
    data = load.colvar(input_col, 'as_pandas')
    data = data[['time', 'pp.proj']].astype({'time': float, 'pp.proj': float})
    # remove any overflow data
    too_long = data[data.time > (max_t*1000)+1].index
    data.drop(too_long, inplace=True)
    # convert projection to Angstroms
    data['pp.proj'] = data['pp.proj'].apply(lambda x: x * 10)
    # calculate status of ligand position: 1 = bound, 2 = unbound
    data['status'] = data['pp.proj'].apply(bound_check, args=([bounds]))
    # remove data without status i.e. not bound or unbound
    middle_ind = data[data.status == 0].index
    data.drop(middle_ind, inplace=True)
    # calculate differences in status column (diff. of 1 = transition)
    data['diffs'] = data.status.diff()
    # identify transitions
    rx_ind = data[data.diffs != 0].index
    # extract times as list
    rx = [int(floor(t)) for t in data.loc[rx_ind[1:]].time.tolist()]
    # count number of recrossings
    N = int(floor(len(rx)/2))

    return N, rx


def dif_multiplot2(pdb_list, cv='proj', rx=False):
    """ Plot 6x2 diffusion plots for paper """
    fig = plt.figure(figsize=(10, 14))
    axes = fig.subplots(6, 2, sharex=True, sharey='row')
    fig.subplots_adjust(left=0.15, hspace=0.2, wspace=0.1)

    A = np.arange(len(pdb_list))
    if '5am1' in pdb_list:
        site = 'tunnel'
    elif '5am4' in pdb_list:
        site = 'RHS'
    elif '5alo' in pdb_list:
        site = 'LHS'
    else:
        print('SITE NOT FOUND')

    if cv == 'both':
        cv_list = ['proj', 'ext']
        n = 0
        for pdb in pdb_list:
            c = 0
            for CV in cv_list:
                f = 0
                for fs in ['FS1', 'FS2']:
                    data = load.colvar('{p}_{f}/{p}-{f}_OLD.colvar'.format(p=pdb, f=fs), 'as_pandas')
                    gr.diffusion(data, pdb, CV, bounds, axes[n+c, f], colours[c2p[pdb]])
                    f += 1
                c += 1
            n += 2
        for i in A:
            fig.text(0.06, 0.225+(A[-i-1]*0.265), pdb_list[i], ha='center', fontsize=14, color=colours[c2p[pdb_list[i]]])
            fig.text(0.1, 0.3+(A[-i-1]*0.265), 'pp.'+cv_list[0]+' / $\mathrm{\AA}$', va='center', rotation='vertical', fontsize=10)
            fig.text(0.1, 0.17+(A[-i-1]*0.265), 'pp.'+cv_list[1]+' / $\mathrm{\AA}$', va='center', rotation='vertical', fontsize=10)

        for i in np.arange(len(cv_list)):
            fig.text(0.33+(i*0.38), 0.9, fun_labels['FS'+str(i+1)], ha='center', fontsize=18, color=colours[i*2])
            fig.text(0.33+(i*0.38), 0.065, 'Simulation Time / ns', ha='center', fontsize=10)
        #fig.legend([str(x)+' $\mathrm{\AA}$' for x in lines] + ['Reweight'], loc='lower center', ncol=6, frameon=False)   
    # plot with recrossings identified and labelled
    elif rx:
        n = 0
        for pdb in pdb_list:
            f = 0
            for fs in ['FS1', 'FS2']:
                load_col = '{p}_{f}/{p}-{f}_OLD.colvar'.format(p=pdb, f=fs)
                data = load.colvar(load_col, 'as_pandas')
                bounds = BOUND_STATES[site][0] if fs=='FS1' else BOUND_STATES[site][1]
                # run recrossing counting function
                N, rx = identify_recross(load_col, bounds)
                # add number of rxs to plot
                axes[n,f].text(0.9, 0.9, 'rx = {}'.format(N), horizontalalignment='center',
                        verticalalignment='center', transform=axes[n,f].transAxes,
                        bbox=dict(facecolor='w', edgecolor='red', alpha=0.8))
                # vertical lines at rx transition points 
                for t in range(N):
                    box = plt.Rectangle((rx[t]/1000 + 5, 5.), (rx[t+1]-rx[t])/1000 -5, 35., ls='--', fc='none', ec='r', alpha=0.7)
                    axes[n,f].add_patch(box)
                    #axes[n,f].axvline(t/1000, ls='--', alpha=0.7, c='r')
                # plot rest of diffusion plot
                gr.diffusion(data, pdb, cv, bounds, axes[n, f], colours[c2p[pdb]])
                f += 1
            n += 1
        # Y axis and plot labels
        for i in A:
            fig.text(0.06, 0.17+(A[-i-1]*0.132), pdb_list[i], ha='center', fontsize=14, color=colours[c2p[pdb_list[i]]])
            fig.text(0.1, 0.17+(A[-i-1]*0.132), 'pp.' + cv + ' / $\mathrm{\AA}$', va='center', rotation='vertical', fontsize=10)
        # X axis and plot labels
        for i in np.arange(2):
            fig.text(0.33+(i*0.38), 0.9, fun_labels['FS'+str(i+1)], ha='center', fontsize=18, color=colours[i*2])
            fig.text(0.33+(i*0.38), 0.065, 'Simulation Time / ns', ha='center', fontsize=10)
    # plot single CV without recrossing info. 
    else:
        n = 0
        for pdb in pdb_list:
            f = 0
            for fs in ['FS1', 'FS2']:
                data = load.colvar('{p}_{f}/{p}-{f}_OLD.colvar'.format(p=pdb, f=fs), 'as_pandas')
                bounds = BOUND_STATES[site][0] if fs=='FS1' else BOUND_STATES[site][1]
                gr.diffusion(data, pdb, cv, bounds, axes[n, f], colours[c2p[pdb]])
                f += 1
            n += 1
        for i in A:
            fig.text(0.06, 0.17+(A[-i-1]*0.132), pdb_list[i], ha='center', fontsize=14, color=colours[c2p[pdb_list[i]]])
            fig.text(0.1, 0.17+(A[-i-1]*0.132), 'pp.' + cv + ' / $\mathrm{\AA}$', va='center', rotation='vertical', fontsize=10)
        for i in np.arange(2):
            fig.text(0.33+(i*0.38), 0.9, fun_labels['FS'+str(i+1)], ha='center', fontsize=18, color=colours[i*2])
            fig.text(0.33+(i*0.38), 0.065, 'Simulation Time / ns', ha='center', fontsize=10)

    fig.savefig('Figures/DIFF_paper_{}.png'.format('_'.join(pdb_list) if not rx else '_'.join(pdb_list)+'_RX'), dpi=300, bbox_inches='tight')

def highlight5am0(cbar_max):
    """ FES over time for 5am0 """
    fig = plt.figure(figsize=(14,14))
    axes = fig.subplots(3,2,sharex=True,sharey=True)
    directory = "./" 
    t_list = [30, 50, 75]

    n = 0
    for t in t_list:
        f = 0
        for fs in ['FS1','FS2']: 
            #plot_num = int(''.join(['2',str(len(pdb_list)),str(n)]))
            #print(plot_num)
            fes_data, axis_labels = load.fes("{}/5am0_{}/fes/fes_{}.dat".format(directory, fs, t), False)   
            cmap = gr.two_cv_contour(fes_data, '5am0', fs, axis_labels, cbar_max, 'OLD_FES', "Figures/{}".format(directory), axes[n,f])
            f+=1
        n += 1

    plt.subplots_adjust(left=0.15)
    
    # Y LABELS
    A = np.arange(len(t_list))
    for i in A:
        fig.text(0.1, 0.24+(A[-i-1]*0.26), 
                 axis_labels[1]+" / $\mathrm{\AA}$", va='center', rotation='vertical', fontsize=10)
        fig.text(0.06, 0.23+(A[-i-1]*0.26), str(t_list[i]*10)+'ns', ha='center', fontsize=14)
    # X LABELS
    for i in [0,1]:
        fig.text(0.34+(i*0.39), 0.065, axis_labels[0]+" / $\mathrm{\AA}$", ha='center', fontsize=10)
        fig.text(0.34+(i*0.39), 0.9, fun_labels['FS'+str(i+1)], ha='center', fontsize=18, color=colours[i*2])
    fig.subplots_adjust(hspace=0.05, wspace=0.05)
    fig.subplots_adjust(right=0.915)
    cax = plt.axes([0.93, 0.11, 0.01, 0.77])
    cbar = plt.colorbar(cmap,cax=cax, aspect=10, ticks=np.arange(0.,cbar_max+1,2.0))
    cbar.set_label("Free Energy / kcal $\mathrm{mol^{-1}}$", fontsize=10) 
    fig.savefig('Figures/5am0_FES_oT.png', dpi=300, bbox_inches='tight')
    
def highlight5am0_dgdt():
    """ dG vs t for 5am0 only"""
    
    fig = plt.figure(figsize=(14,6))
    axes = fig.subplots(1, 2, sharey=True)
    errors = pd.read_csv('5am0_extended.csv', sep=',')
    #DATA_FILE = "./dG_database.h5"
    #database = pd.read_hdf(DATA_FILE, key='deltag') 
    for i in [0,1]:
        data_file = open("5am0_FS{}_dGval_750ns.p".format(i+1),"rb")
        y_data = pickle.load(data_file)
        print('FS{}'.format(i+1))
        for val in np.arange(len(y_data)):
            print('{:4d}  {:6.3f}'.format(val, y_data[val]))
    
        gr.dgdt(y_data, -8.2, axes[i], err=errors['FS{}'.format(i+1)])
        fig.text(0.32+(i*0.4), 0.01, "Simulation Time / ns", ha='center', fontsize=10)
        fig.text(0.32+(i*0.4), 0.9, fun_labels['FS'+str(i+1)], ha='center', fontsize=18, color=colours[i*2])
        fig.text(0.08, 0.5, '$\Delta$G / kcal $\mathrm{mol^{-1}}$', va='center', rotation='vertical', fontsize=10)
    fig.subplots_adjust(hspace=0.05, wspace=0.05)
    fig.subplots_adjust(right=0.915)
    fig.savefig('Figures/5am0_dG_oT.png', dpi=300, bbox_inches='tight')

def dgdt_plots(pdb_list):
    """ dG vs t """

    fig = plt.figure(figsize=(10, 14))
    axes = fig.subplots(len(pdb_list), 2, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.15, hspace=0.2, wspace=0.1)

    DATA_FILE = "./dG_database.h5"
    database = pd.read_hdf(DATA_FILE, key='deltag')
    exp_vals = pd.read_csv('dg_values.csv', sep=',')
    rhs_errs = pd.read_csv('ERR_oT_Frhs.csv', sep=',')
    lhs_errs = pd.read_csv('ERR_oT_Flhs.csv', sep=',')

    A = np.arange(len(pdb_list))
    for j in A:
        y_data = [database[pdb_list[j]].FS1.to_numpy(), 
                  database[pdb_list[j]].FS2.to_numpy()]
        print(min(y_data[0]), min(y_data[1]))

        dg_exp = float(exp_vals[exp_vals['pdb'] == pdb_list[j]].exp)

        errors = [rhs_errs[pdb_list[j]], lhs_errs[pdb_list[j]]]

        print(pdb_list[j], dg_exp)
        # Y axis and plot labels
        fig.text(0.06, 0.17+(A[-j-1]*0.132), pdb_list[j], ha='center', fontsize=14, color=colours[c2p[pdb_list[j]]])
        fig.text(0.1, 0.17+(A[-j-1]*0.132), '$\Delta$G / kcal $\mathrm{mol^{-1}}$', va='center', rotation='vertical', fontsize=10)
     
        # X axis and plot labels
        for i in np.arange(2):
            gr.dgdt(y_data[i], dg_exp, axes[j, i], err=errors[i])
            fig.text(0.33+(i*0.38), 0.9, fun_labels['FS'+str(i+1)], ha='center', fontsize=18, color=colours[i*2])
            fig.text(0.33+(i*0.38), 0.065, 'Simulation Time / ns', ha='center', fontsize=10)

        #fig.text(0.1, 0.24+(A[-j-1]*0.26), 
                 #'$\Delta$G / kcal $\mathrm{mol^{-1}}$', va='center', rotation='vertical', fontsize=10)
        #fig.text(0.06, 0.23+(A[-j-1]*0.26), pdb_list[j], ha='center', fontsize=14)
 
    fig.savefig('Figures/dG_oT_'+'_'.join(pdb_list)+'.png', dpi=300, bbox_inches='tight')

def ddGdt_plot(n_reps, pdb_clust, ts, highlight=None):
    """ ddG vs t """ 
    #colours = ['#7e57c2','#3f51b5','#1e88e5','#0277bd','#00838f','#004d40']  # dark green to purple
    #colours = ['#061e46','#0d47a1','#1976d2','#2196f3','#64b5f6','#bbdefb']  # blues
    colours = ['xkcd:red','xkcd:orange', 'xkcd:yellow', 'xkcd:green', 'xkcd:blue', 'xkcd:violet']


    if len(ts) == 1:
        fig, ax = plt.subplots(figsize=(5.5,5))
    else:
        fig, ax = plt.subplots(figsize=(9.6,5))
        ax.axvline(x=300, ymin=-5., ymax=12.,  color='#0d47a1', ls='--')
    exp_vals = pd.read_csv('dg_values.csv', sep=',')
    
    first = True 
    for T in ts:
        data = {}
        over = np.arange(n_reps) if highlight is None else highlight
        for r in over:
            csv = '../../SWISH/ANALYSIS/SWISH_dGoT/{}ns/SWISH_dGoT_rep{}.csv'.format(T,r)
            database = pd.read_csv(csv, sep=',')
            database = database.filter(pdb_clust)
            for c in database.columns:
                database[c] = database[c] - float(exp_vals[exp_vals['pdb']==c].exp)
            database = database.abs()
            data[r] = [database.mean(axis=1), database.std(axis=1)]
            
            max_ylen = int(T/10)
            x = np.linspace(0., max_ylen*10, max_ylen+1)
            if first: 
                ax.plot(x, data[r][0], c=colours[r], zorder=1, label='$\mathrm{\lambda_{%s}}$' % r)
                if highlight is not None:
                    ax.fill_between(x, data[r][0]+data[r][1], data[r][0]-data[r][1],facecolor=colours[r], alpha=.1)
            else:
                upto = int((T - ts[ts.index(T)-1])/10)
                ax.plot(x[-upto:], data[r][0][-upto:], c = colours[r], zorder=1,)
                if highlight is not None:
                    ax.fill_between(x[-upto:], data[r][0][-upto:]+data[r][1][-upto:], data[r][0][-upto:]-data[r][1][-upto:],facecolor=colours[r], alpha=.1)
        first=False
        print(data)
    
    ax.axhline(y=0., xmin=0., xmax=max(x), c='xkcd:green', ls='--') 
    ax.axhline(y=2.0, xmin=0., xmax=max(x), c='xkcd:green', alpha=0.2) 
    ax.axhspan(ymin=0., ymax=2., xmin=0., xmax=max(x),  facecolor='xkcd:green', alpha=0.1) 
    ax.axhline(y=3.5, xmin=0., xmax=max(x),  color='xkcd:orange', alpha=0.2)
    ax.axhspan(ymin=0., ymax=3.5, xmin=0., xmax=max(x),  facecolor='xkcd:orange', alpha=0.1)
    ax.legend()
    

    ax.set_xlim([-5., max(ts)])
    if highlight is None:
        ax.set_ylim([-0.5,11.])
    else: 
        ax.set_ylim([-2,20.])
    
    ax.set_xlabel("Simulation Time / ns", ha='center', fontsize=10) 
    ax.set_ylabel('$\Delta\Delta$G / kcal $\mathrm{mol^{-1}}$', va='center', rotation='vertical', fontsize=10) 
          
    fig.savefig('Figures/ddG_v_T_{}_{}{}.png'.format(pdb_clust[0], max(ts), '' if highlight is None else '_H'), dpi=300, bbox_inches='tight') 

def SWISH_diff_multi(pdb_list, size, demux=False, rx=False, boxes=False):
    """ Plot 6x2 diffusion plots for paper """  

    _lambda = 1

    if '5am3' in pdb_list:
        site = 'tunnel'
    elif '5aia' in pdb_list:
        site = 'RHS'
    elif '5akk' in pdb_list:
        site = 'LHS'
    else:
        print('SITE NOT FOUND')
    
    fig = plt.figure(figsize=(size[0]*4,size[1]*4.5) if size == [3, 2] else (12, 6.5))
    axes = fig.subplots(size[0],size[1],sharex=False,sharey=True) 
    axes = axes.ravel()
    A = np.arange(len(pdb_list)) 
    micros = ['5alg', '5alh', '5aly', '5am0', '5am3']
    for n in A:
        pdb = pdb_list[n]
        #colvar_path = "/media/arc/cem/UCB_PROJECT/SWISH/ANALYSIS/300ns/{p}_FS2/colvars/COLVAR.{l}".format(p=pdb, l=_lambda)\
            #if not demux else "/media/arc/cem/UCB_PROJECT/SWISH/ANALYSIS/300ns/{p}_FS2/demux/colvars/{l}_{p}_300ns_demux.colvar".format(l=_lambda, p=pdb)
        if pdb in micros:
            T = 1000
        elif pdb == "5alo":
            T = 500
        else:
            T = 300

        if not demux:
            colvar_path = "/media/rhys/storage/UCB_Project/SWISH/ANALYSIS/{t}ns/{p}_FS2/colvars/COLVAR.{l}".format(p=pdb, l=_lambda, t=T)
        else:
            colvar_path = "/media/rhys/storage/UCB_Project/SWISH/ANALYSIS/{t}ns/{p}_FS2/demux/colvars/{l}_{p}_{t}ns_demux.colvar".format(p=pdb, l=_lambda, t=T)
        colvar_path = "./COMetPath_Colvars/{p}/COLVAR_{p}_short".format(p=pdb)

        print(colvar_path)
        data = load.colvar(colvar_path, 'as_pandas')

        bounds = BOUND_STATES[site][1]
        if rx:
            # run recrossing counting function
            N, rx = identify_recross(colvar_path, bounds)
            # add number of rxs to plot
            axes[n].text(0.89, 0.9, 'rx = {}'.format(N), horizontalalignment='center',
                    verticalalignment='center', transform=axes[n].transAxes,
                    bbox=dict(facecolor='w', edgecolor='red', alpha=0.8))
            if boxes:
                for t in range(N):
                    box = plt.Rectangle((rx[t]/1000 + 5, 5.), (rx[t+1]-rx[t])/1000 -5, 35., ls='--', fc='none', ec='r', alpha=0.7)
                    axes[n].add_patch(box)

        lin_col = colours[c2p[pdb]]
        label_N = (size[0]*size[1]) - size[1] - 1
        gr.SWISH_diffusion(pdb, _lambda, axes[n], data, lin_col, n, label_N, bounds, demux)
    k = 0
    for j in np.arange(size[0]):
        for i in np.arange(size[1]):
            lin_col = colours[c2p[pdb_list[k]]]
            if size == [3, 2]:
                #        X disp             Y disp
                fig.text(0.135+(i*0.42), 0.3+(np.arange(size[0])[-j-1]*0.272), pdb_list[k], ha='left', fontsize=20, color=lin_col)
            if size == [2, 2]:
                fig.text(0.13+(i*0.4), 0.43+(np.arange(size[0])[-j-1]*0.395), pdb_list[k], ha='left', fontsize=20, color=lin_col)
            k += 1
    #fig.subplots_adjust(hspace=0.05, wspace=0.05) 
    s = pdb
    if demux:
        s += '_demux'
    if rx:
        s += '_RX'
    fig.savefig('Figures/COMET_diffMULTI_{}.png'.format(s), dpi=300, bbox_inches='tight')
    plt.close(fig)

def SWISH_energies(pdb, n_res):
    e_path = "/media/arc/cem/UCB_PROJECT/SWISH/ANALYSIS/300ns/5akk_FS2/energies/"
    binwidth = 100
    fig = plt.figure(figsize=(8, 6))
    ax =  fig.subplots(1,1)

    for i in np.arange(n_res):
        energy_data = load.xvg("{}{}_{}_300ns_energy.xvg".format(e_path, i, pdb))
        data = energy_data[1]
        ax.hist(data, bins=np.arange(min(data), max(data) + binwidth, binwidth))

    
    fig.savefig('Figures/SWISH_energies_{}.png'.format(pdb), dpi=300, bbox_inches='tight')

def SWISH_exchange(pdb, n_res):
    sns.set(color_codes=True)
    #o, 'as_pandas')ut_dir = 'Figures/{}_exchanges/'.format(pdb)
    if not os.path.exists(out_dir): 
        os.mkdir(out_dir) 
    index_file = '/media/arc/cem/UCB_PROJECT/SWISH/C-TERM/{}/SCALE_5/replica_index.xvg'.format(pdb)
    df=pd.read_csv(index_file, sep='\s+', header=None, names=['time', 'r0', 'r1', 'r2', 'r3','r4', 'r5','r6', 'r7'])
    df.time = df.time/1000

    for i in np.arange(n_res):
        print("INFO:    Starting replica {}".format(i))
        g=sns.jointplot(x="time", y="r{}".format(i), data=df, kind="kde", color="m")
        g.savefig('{}{}_r{}.png'.format(out_dir, pdb, i), transparent=True, dpi=300, bbox_inches='tight') 

def SWISH_fes(pdb, n_res, time, cbar_max):
    fig = plt.figure(figsize=(18,10))
    axes = fig.subplots(2,3,sharex=True,sharey=True)
    axes = axes.ravel() 
    directory = "/media/arc/cem/UCB_PROJECT/SWISH/ANALYSIS/{t}ns/{p}_FS2/FES/".format(t=time, p=pdb)  
    axis_labels = ['pp.proj', 'pp.ext'] #if not is_rew else ['RMSD OUT', 'RMSD IN']
    for REP in np.arange(n_res):
        #plot_num = int(''.join(['2',str(len(pdb_list)),str(n)]))
        #print(plot_num)
        fes_data = load.fes_simple("{d}{n}_{p}_FS2_{t}ns.fes".format(d=directory, p=pdb, n=REP, t=time), False)   
        cmap = gr.two_cv_contour(fes_data, pdb, 2, axis_labels, cbar_max, 'SWISH_FES', "Figures/{}".format(directory), axes[REP])
        #axes[REP].set_title('Rep. {}'.format(REP))

    plt.subplots_adjust(left=0.15)
    # Y LABELS
    A = np.arange(2)
    n = 0
        
        #fig.text(0.06, 0.23+(A[-i-1]*0.26), pdb_list[i], ha='center', fontsize=14)
    for j in A:     
        fig.text(0.1, 0.34+(A[-j-1]*0.39), axis_labels[1]+" / $\mathrm{\AA}$", va='center', rotation='vertical', fontsize=12)
        for i in [0,1, 2]:
            #        X disp             Y disp
            fig.text(0.375+(i*0.259), 0.46+(A[-j-1]*0.405), '$\mathrm{\lambda_{%s}}$' % str(n), ha='center', fontsize=20)
            n += 1
        
            fig.text(0.28+(i*0.259), 0.065, axis_labels[0]+" / $\mathrm{\AA}$", ha='center', fontsize=12)

        #fig.text(0.34+(j*0.39), 0.9, fun_labels['FS'+str(j+1)], ha='center', fontsize=18, color=colours[j*2])
    
    fig.text(0.08, 0.5, pdb, ha='center', fontsize=18, color=colours[c2p[pdb]])
    fig.subplots_adjust(hspace=0.05, wspace=0.05, right=0.915, top=0.9)
    cax = plt.axes([0.93, 0.11, 0.01, 0.77])
    cbar = plt.colorbar(cmap,cax=cax, aspect=10, ticks=np.arange(0.,cbar_max+1,2.0))
    cbar.set_label("Free Energy / kcal $\mathrm{mol^{-1}}$", fontsize=12) 
    fig.savefig('Figures/SWISH_FES_{}_{}ns.png'.format(pdb, time), dpi=300, bbox_inches='tight')

def autolabel(rects, axis, horizontal=False, z=20):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        if horizontal:
            width = rect.get_width()
            axis.text(rect.get_x() + rect.get_width() / 2, rect.get_y() + rect.get_height()/2.,
                    '%.2f' % width,
                    ha='center', va='center', color='white', zorder=z)
        else:
            height = rect.get_height()
            axis.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom', zorder=z) 

def performance_bar():
 
    # change colours 
    
    labels = ['4x Platinum 8160 (96 cores)', 
                '8x Platinum 8160 (192 cores)', 
                '12x Platinum 8160 (288 cores)', 
                '1x Tesla P100',
                '1x Tesla P100',
                '6x Tesla P100',
             ]

    pos = [(1, 5), 
           (2, 4),
           (0, 3)
           ]
    vals = [52.5,
            60.7,
            60.0,
            42.0,
            54.0,
            72.1,
            ]
    scale = 0.3
    x = np.arange(len(labels))*scale
    width = 0.24
    
    fig, ax = plt.subplots(figsize=(9,4)) 

    bar1 = ax.barh([y*scale for y in pos[0]], [vals[i] for i in pos[0]], width, label='Fun-metaD', color='#5a7d9a') 
    bar2 = ax.barh([y*scale for y in pos[1]], [vals[i] for i in pos[1]], width, label='Fun-SWISH', color='#0a888a') 
    bar3 = ax.barh([y*scale for y in pos[2]], [vals[i] for i in pos[2]], width, label='COMetPath', color='#1e488f')  
    ax.axhline(y=2.5*scale, xmin=0., xmax=max(vals)+10, c='xkcd:charcoal', ls='--')
   
    ax.set_xlim([0., 80.])
    ax.set_xlabel('Simulation Performance / ns·$\mathrm{day^{-1}}$', fontsize=9)
    ax.set_yticks(x)
    ax.set_yticklabels(labels, fontsize=8)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=8) 
    
    fig.subplots_adjust(bottom=0.2) 
    ax.legend( loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, 0), bbox_transform=fig.transFigure)
    
    autolabel(bar1, ax, True)
    autolabel(bar2, ax, True)
    autolabel(bar3, ax, True)
    hpc = ['CPU', 'GPU']
    A = np.arange(2)
    for i in A:
        fig.text(0.88, 0.48+(A[-i-1]*0.34), hpc[i], ha='right', fontsize=14, color='xkcd:charcoal')

    plt.gca().invert_yaxis()

    fig.savefig('Figures/Performance.png', dpi=300, bbox_inches="tight")

def performance_bar2():
 
    # change colours 
    
    labels = ['Fun-metaD',
              'Fun-SWISH',
              'COMetPath',
              ] * 2

 
    vals = [np.array([35.56, 32.9, 25.1, 38.5, 54.0, 42.0]),
            #np.array([50.06, 46.48, 50.8, 59.8, 0.0, 0.0]),
            np.array([50.06, 46.48, 50.8, 0.0, 0.0, 0.0]),
            np.array([53.3, 51.6, 62.4, 59.8, 0.0, 0.0]),
            np.array([61.6, 53.9, 0.0, 70.1, 0.0, 0.0]),
            np.array([64.8, 57.4, 0.0, 0.0, 0.0, 0.0]),
            np.array([66.2, 60.4, 0.0, 72.9, 0.0, 0.0]),
            ]
 
    print(labels)
    ind = np.arange(len(vals[0]))
    print(ind)

    fig, ax = plt.subplots(figsize=(9,4))

    BARS = [None]*6
    BARS[5] = ax.barh(ind, vals[5]-vals[4], label='6 Nodes', alpha=0.25, left=vals[4], zorder=1)
    BARS[4] = ax.barh(ind, vals[4]-vals[3], label='5 Nodes', alpha=0.40, left=vals[3], zorder=2)
    BARS[3] = ax.barh(ind, vals[3]-vals[2], label='4 Nodes', alpha=0.55, left=vals[2], zorder=3)
    BARS[2] = ax.barh(ind, vals[2]-vals[1], label='3 Nodes', alpha=0.70, left=vals[1], zorder=4)
    BARS[1] = ax.barh(ind, vals[1]-vals[0], label='2 Nodes', alpha=0.85, left=vals[0], zorder=5)
    BARS[0] = ax.barh(ind, vals[0], label='1 Node', alpha=1.0, zorder=10)

    bar_colours = ['#5a7d9a', '#0a888a', '#1e488f'] * 2

    for bar in BARS:
        for i in np.arange(6):
            bar[i].set_color(bar_colours[i])

    ax.axhline(y=2.5, xmin=0., xmax=85, c='xkcd:charcoal', ls='--')
   
    ax.set_xlim([0., 85])
    ax.set_xlabel('Simulation Performance / ns·$\mathrm{day^{-1}}$', fontsize=9)
    x = np.arange(len(labels))
    ax.set_yticks(x)
    ax.set_yticklabels(labels, fontsize=8)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=8)
    
    fig.subplots_adjust(bottom=0.2)
    H, L = ax.get_legend_handles_labels()
    ax.legend(H[::-1], L[::-1], loc='lower center', ncol=6, frameon=False, bbox_to_anchor=(0.5, 0), bbox_transform=fig.transFigure)
    
    autolabel(BARS[0], ax, True)
    hpc = ['CPU', 'GPU']
    A = np.arange(2)
    for i in A:
        fig.text(0.88, 0.48+(A[-i-1]*0.34), hpc[i], ha='right', fontsize=14, color='xkcd:charcoal')

    plt.gca().invert_yaxis()

    fig.savefig('Figures/Performance.png', dpi=300, bbox_inches="tight")

###############################################################################

if __name__ == "__main__":
    
    #### FES MULTIPLOTS ####
    pdb_orders = [['5am1', '5am3', '5alg'],
                  ['5alp', '5alh', '5ai5'],
                  ['5am4', '5aly', '5akh'],
                  ['5am0', '5alx', '5aia'],
                  ['5alo', '5alt', '5akg'],
                  ['5akk', '5ai0', '5ak6'], 
                  ]
    core_list = ['5alp', '5alg',
                 '5alx', '5aia',
                 '5akk', '5alt',
                 ]
    
    #for pdb_list in pdb_orders: 
        #fes_multiplot(pdb_list, cbar_max=22)

    #### FES HIGHLIGHTS ####
    #fes_highlight('5ak6')


    
    #### CONVERGENCE MULTIPLOTS ####
    #for i in [0,2,4]:
        #pdb_cut = pdb_orders[i].copy()
        #pdb_cut.extend(pdb_orders[i+1])    
        #conv_multiplot(pdb_cut, 'cut')
   #  conv_multiplot(core_list,'cut')

    #### DIFFUSION MULTIPLOTS ####
    #dif_multiplot2(['5aly','5alp','5alt'], cv='both')
    
    # dif_multiplot2(core_list, cv='proj')

    #### 5AM0 PLOTS ####
    #highlight5am0(cbar_max=22)


    #### dG vs. TIME PLOTS ####
    #for l in pdb_orders:
        #dgdt_plots(l)
        #dif_multiplot2(l)

    
    #for i in [0,2,4]:
        #pdb_cut = pdb_orders[i].copy()
        #pdb_cut.extend(pdb_orders[i+1])  
        #ddGdt_plot(6, pdb_cut, [300]) 
    #ddGdt_plot(6, ['5am1', '5am3', '5alg', '5alp', '5alh', '5ai5'], [300, 1000])
    #ddGdt_plot(6, ['5am1', '5am3', '5alg', '5alp', '5alh', '5ai5'], [300, 500], highlight=[1,5]) 
    

    #### SWISH PLOTS ####
    #plist = ['5akk', '5aly','5alp', '5am1', '5am4', '5alt']
    #for plist in pdb_orders:
    #for pdb in plist:
        #SWISH_diffusion(pdb, 1)
        #SWISH_diffusion(pdb, 1, demux=True)
        #SWISH_fes(pdb, 6, 300, 22)
    #SWISH_fes('5aia', 6, 300, 22)
    #SWISH_fes('5am3', 6, 500, 22)
    #SWISH_fes('5am3', 6, 1000, 22)

    #SWISH_diff_multi(['5alp', '5alg', '5alx', '5aia', '5akk', '5alt'], [3, 2])
    #SWISH_diff_multi(['5am3', '5alh', '5aly', '5am0'], [2, 2])

    #SWISH_energies('5akk', 6)
    #SWISH_exchange('5akk',6) 

    for i in [0,2,4]:
        pdb_cut = pdb_orders[i].copy()
        pdb_cut.extend(pdb_orders[i+1])
        #dif_multiplot2(pdb_cut, cv='proj')
        #dif_multiplot2(pdb_cut, cv='proj', rx=True)
        #SWISH_diff_multi(pdb_cut, [3, 2])
        #SWISH_diff_multi(pdb_cut, [3, 2], rx=True)

        #dgdt_plots(pdb_cut)
    #highlight5am0_dgdt()
    #SWISH_diff_multi(['5am3','5alg','5alp','5am3','5alg','5am3'], [3,2], rx=True,boxes=True)
    #SWISH_diff_multi(['5aia','5alx','5aly','5aia','5alx','5aly'], [3,2], rx=True,boxes=True)
    #SWISH_diff_multi(['5akk','5akg','5alt','5akk','5akg','5alt'], [3,2], rx=True,boxes=True)
    performance_bar2()
