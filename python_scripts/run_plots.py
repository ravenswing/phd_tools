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
            fes_data, axis_labels = load.fes("{}/{}-{}_{}.fes".format(directory, pdb, fs, fes_type), is_rew)
            cmap = gr.two_cv_contour(fes_data, pdb, fs, axis_labels, cbar_max, 'OLD_{}'.format(fes_type), "Figures/{}".format(directory), axes[f,n])
            f += 1
        n += 1
    plt.subplots_adjust(left=0.15)
    labels = {'OLD':['Raw Data', 'pp.proj / $\mathrm{\AA}$'],
              'REW':['Reweighted', axis_labels[0]+' / $\mathrm{\AA}$'],
             }
    # Y LABELS
    for i in [0, 1]:
        fig.text(0.04, 0.3+(i*0.38), 'FS'+str(i+1), ha='center', fontsize=14)
        fig.text(0.08, 0.31+(i*0.38), 'pp.ext / $\mathrm{\AA}$', va='center', rotation='vertical', fontsize=10)

    # X LABELS
    for i in np.arange(len(fes_list)):
        fig.text(0.33+(i*0.41), 0.9, labels[fes_list[i]][0], ha='center', fontsize=14)
        fig.text(0.33+(i*0.41), 0.065, labels[fes_list[i]][1], ha='center', fontsize=10)

    fig.subplots_adjust(hspace=0.05, wspace=0.15)
    fig.subplots_adjust(right=0.915)
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

def dif_multiplot2(pdb_list, cv='proj'):
    """ Plot 6x2 diffusion plots for paper """
    fig = plt.figure(figsize=(10, 14))
    axes = fig.subplots(6, 2, sharex=True, sharey='row')
    fig.subplots_adjust(left=0.15, hspace=0.2, wspace=0.1)

    A = np.arange(len(pdb_list))
    
    if cv == 'both':
        cv_list = ['proj', 'ext']
        n = 0
        for pdb in pdb_list:
            c = 0
            for CV in cv_list:
                f = 0
                for fs in ['FS1', 'FS2']:
                    data = load.colvar('{p}_{f}/{p}-{f}_OLD.colvar'.format(p=pdb, f=fs), 'as_pandas')
                    gr.diffusion(data, CV, axes[n+c, f], colours[c2p[pdb]])
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
    else:
        n = 0
        for pdb in pdb_list:
            f = 0
            for fs in ['FS1', 'FS2']:
                data = load.colvar('{p}_{f}/{p}-{f}_OLD.colvar'.format(p=pdb, f=fs), 'as_pandas')
                gr.diffusion(data, cv, axes[n, f], colours[c2p[pdb]])
                f += 1
            n += 1
        for i in A:
            fig.text(0.06, 0.17+(A[-i-1]*0.132), pdb_list[i], ha='center', fontsize=14, color=colours[c2p[pdb_list[i]]])
            fig.text(0.1, 0.17+(A[-i-1]*0.132), 'pp.' + cv + ' / $\mathrm{\AA}$', va='center', rotation='vertical', fontsize=10)
        for i in np.arange(2):
            fig.text(0.33+(i*0.38), 0.9, fun_labels['FS'+str(i+1)], ha='center', fontsize=18, color=colours[i*2])
            fig.text(0.33+(i*0.38), 0.065, 'Simulation Time / ns', ha='center', fontsize=10)

    fig.savefig('Figures/DIFF_paper_'+'_'.join(pdb_list)+'.png', dpi=300, bbox_inches='tight')

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
    #DATA_FILE = "./dG_database.h5"
    #database = pd.read_hdf(DATA_FILE, key='deltag') 
    for i in [0,1]:
        data_file = open("5am0_FS{}_dGval_750ns.p".format(i+1),"rb")
        y_data = pickle.load(data_file)
        print('FS{}'.format(i+1))
        for val in np.arange(len(y_data)):
            print('{:4d}  {:6.3f}'.format(val, y_data[val]))
    
        gr.dgdt(y_data, -8.2, axes[i])  
        fig.text(0.32+(i*0.4), 0.01, "Simulation Time / ns", ha='center', fontsize=10)
        fig.text(0.32+(i*0.4), 0.9, fun_labels['FS'+str(i+1)], ha='center', fontsize=18, color=colours[i*2])
        fig.text(0.08, 0.5, '$\Delta$G / kcal $\mathrm{mol^{-1}}$', va='center', rotation='vertical', fontsize=10)
    fig.subplots_adjust(hspace=0.05, wspace=0.05)
    fig.subplots_adjust(right=0.915)
    fig.savefig('Figures/5am0_dG_oT.png', dpi=300, bbox_inches='tight')

def dgdt_plots(pdb_list):
    """ dG vs t for 5am0 only"""
    
    fig = plt.figure(figsize=(14,14))
    axes = fig.subplots(3,2,sharex=True,sharey=True)
     
    DATA_FILE = "./dG_database.h5"
    database = pd.read_hdf(DATA_FILE, key='deltag')
    exp_vals = pd.read_csv('dg_values.csv', sep=',')

    A = np.arange(len(pdb_list)) 
    for j in A:
        y_data = [database[pdb_list[j]].FS1.get_values(), 
                  database[pdb_list[j]].FS2.get_values()]
        #print(min(y_data[0]), min(y_data[1]))
            
        dg_exp = float(exp_vals[exp_vals['pdb']==pdb_list[j]].exp)
        
        print(pdb_list[j], dg_exp)

        for i in [0,1]:
            gr.dgdt(y_data[i], dg_exp, axes[j, i])    
            fig.text(0.34+(i*0.39), 0.065, "Simulation Time / ns", ha='center', fontsize=10)
            fig.text(0.34+(i*0.39), 0.9, fun_labels['FS'+str(i+1)], ha='center', fontsize=18, color=colours[i*2])

        fig.text(0.1, 0.24+(A[-j-1]*0.26), 
                 '$\Delta$G / kcal $\mathrm{mol^{-1}}$', va='center', rotation='vertical', fontsize=10)
        fig.text(0.06, 0.23+(A[-j-1]*0.26), pdb_list[j], ha='center', fontsize=14)
         
    fig.subplots_adjust(left=0.15, right=0.915, hspace=0.08, wspace=0.05) 
    fig.savefig('Figures/dG_oT_'+'_'.join(pdb_list)+'.png', dpi=300, bbox_inches='tight')

def SWISH_diffusion(pdb, _lambda, demux=False):
    """ Plot 6x2 diffusion plots for paper """ 
    c2p = {'5aly':0, '5akk':2, '5alp':1, '5ai0':2, '5ai5':1, '5alt':2, '5am4':0, '5am1':1}

    fig = plt.figure(figsize=(8, 6))
    ax =  fig.subplots(1,1)
    colvar_path = "/media/arc/cem/UCB_PROJECT/SWISH/ANALYSIS/{p}_FS2/colvars/COLVAR.{l}".format(p=pdb, l=_lambda)\
        if not demux else "/media/arc/cem/UCB_PROJECT/SWISH/ANALYSIS/{p}_FS2/demux/colvars/{l}_{p}_300ns_demux.colvar".format(l=_lambda, p=pdb)
    data = load.colvar(colvar_path, 'as_pandas')
    
    N = 50 if not demux else 1
    data = data.iloc[::N, :]
    cv='proj'
    lin_col=colours[c2p[pdb]]
 
    data['pp.'+cv] = data['pp.'+cv].astype(float)
    data['mean'] = data['pp.'+cv].rolling(25,center=True).mean()
    y1 = data['pp.'+cv].values*10
    y2 = data['mean'].values*10
    x = np.arange(len(y1))*(0.002*N) if not demux else np.arange(len(y1))
    ax.plot(x, y1, c=lin_col, alpha=0.3, lw=0.5)
    ax.plot(x, y2, c='k', alpha=1., lw=1.)
   
    ax.set_ylim([-3, 50.])
    ax.set_yticks(np.arange(6)*10)
    ax.axhline(y=30, xmax=500., c='k', alpha=0.5, lw=1., ls='--')
    ax.axhline(y=7, xmax=500., c='k', alpha=0.5, lw=1., ls='--')
     
    ax.set_xlim([0, 300.])
    ax.set_xticks(np.linspace(0., 300., num=7))
    ax.grid(alpha=0.3)
    
    ax.set_ylabel('pp.proj / $\mathrm{\AA}$', va='center', rotation='vertical', labelpad=10, fontsize=10)
    ax.set_xlabel('Simulation Time / ns', ha='center', fontsize=10, labelpad=10) 
    
    fig.savefig('Figures/SWISH_dif_{}.png'.format(pdb+'_demux' if demux else pdb), dpi=300, bbox_inches='tight')

def SWISH_energies(pdb, n_res):
    e_path = "/media/arc/cem/UCB_PROJECT/SWISH/ANALYSIS/5akk_FS2/energies/"
    fig = plt.figure(figsize=(8, 6))
    ax =  fig.subplots(1,1)

    for i in np.arange(n_res):
        energy_data = load.xvg("{}{}_{}_300ns_energy.xvg".format(e_path, i, pdb))
        ax.hist(energy_data[1], 50)

    
    fig.savefig('Figures/SWISH_energies_{}.png'.format(pdb), dpi=300, bbox_inches='tight')

def SWISH_exchange(pdb, n_res):
    sns.set(color_codes=True)
    out_dir = 'Figures/{}_exchanges/'.format(pdb)
    if not os.path.exists(out_dir): 
        os.mkdir(out_dir) 
    index_file = '/media/arc/cem/UCB_PROJECT/SWISH/C-TERM/{}/SCALE_5/replica_index.xvg'.format(pdb)
    df=pd.read_csv(index_file, sep='\s+', header=None, names=['time', 'r0', 'r1', 'r2', 'r3','r4', 'r5','r6', 'r7'])
    df.time = df.time/1000

    for i in np.arange(n_res):
        print("INFO:    Starting replica {}".format(i))
        g=sns.jointplot(x="time", y="r{}".format(i), data=df, kind="kde", color="m")
        g.savefig('{}{}_r{}.png'.format(out_dir, pdb, i), transparent=True, dpi=300, bbox_inches='tight') 


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
    
    for pdb_list in pdb_orders:
        print("")
        #fes_multiplot(pdb_list, cbar_max=22)


    
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
    highlight5am0(cbar_max=22)

    #highlight5am0_dgdt()

    #### dG vs. TIME PLOTS ####
    #for l in pdb_orders:
        #dgdt_plots(l)
        #dif_multiplot2(l)
    
   # dgdt_plots(['5alp', '5alx', '5akk'])

    #### SWISH PLOTS ####
    #plist = ['5akk', '5aly','5alp', '5am1', '5am4', '5alt']
    #for pdb in plist:
        #SWISH_diffusion(pdb, 1)
        #SWISH_diffusion(pdb, 1, demux=True)
    #SWISH_energies('5akk', 6)
    #SWISH_exchange('5akk',6) 
    
    #### FES HIGHLIGHTS ####
    #fes_highlight('5akk')
