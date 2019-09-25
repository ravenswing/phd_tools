"""
===============================================================================
                                    PLOTTING RUN

===============================================================================
"""

import glob
import subprocess
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# local files
import graphics as gr
import load_data as load

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

    for directory in dir_list:
        # extract system info.
        pdb, fs = directory.split("_")
        # make directory for plots
        subprocess.call("mkdir Figures/{}/".format(directory), shell=True)
'''
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
'''

#all_plots(data_dirs)
pdb_list = ['5akk','5alt', '5ai0']
def fes_multiplot(pdb_list):

    fig = plt.figure(figsize=(8,10))
    axes = fig.subplots(3,2,sharex=True,sharey=True)
    n = 0
    for pdb in pdb_list:
        f = 0
        for fs in ['FS1','FS2']:
            directory = '{}_{}'.format(pdb, fs)
            plot_num = int(''.join(['2',str(len(pdb_list)),str(n)]))
            print(plot_num)

            fes_data, axis_labels = load.fes("{}/{}-{}_OLD.fes".format(directory, pdb, fs), False)  
            #fes_data, axis_labels = load.fes("{}/{}-{}_OLD.fes".format(directopdb,fs), False)
            gr.two_cv_contour(fes_data, pdb, fs, axis_labels, 31, 'OLD_FES', "Figures/{}".format(directory), axes[n,f])
            f+=1

        n += 1


    plt.subplots_adjust(left=0.2)
    for i in np.arange(len(pdb_list)):
        fig.text(0.1, 0.165+(i*0.33), axis_labels[1]+' / nm', va='center', rotation='vertical')
        fig.text(0.04,0.165+(i*0.33), pdb_list[i], ha='center')

    fig.text(0.5, 0.04, axis_labels[0]+' / nm', ha='center')
    fig.text(0.33, 0.90, 'FS1', ha='center')
    fig.text(0.66, 0.90, 'FS2', ha='center')
    fig.subplots_adjust(hspace=0.1, wspace=0.1)
    plt.show()
    fig.savefig('Figures/TEST_MULTIPLOT.png', dpi=300, bbox_inches='tight')


fes_multiplot(pdb_list)
