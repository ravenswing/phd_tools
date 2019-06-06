import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy as cp
import pickle
import os
import glob
import sys
import seaborn as sns

###############################################################################
#                                   INPUTS
###############################################################################

# Source directory for the files on your local system
srcdir = '/media/rhys/ExtHD/Project/carlos_peptides/LONG/hydrophilic_brush/helical_brush/'

mol = 'EK5_helix'
figure_path = '{}figures/{}/'.format(srcdir,mol)
pep_length=10
num_pep=9
# default_x = np.linspace(0,2500000.0,num=25001)
plt.style.use('ggplot')

###############################################################################
#                               HEATMAP PLOT FUNC
###############################################################################

def plotHeatMap(systems, limit, distances, coordinates, marker_size, xlabel, ylabel, colourlabel, saveFig, save_folder):
    iad = systems
    cm = plt.cm.get_cmap('plasma')
    m, axis = plt.subplots(2, 2, figsize=(17,14))
    axis = axis.ravel()
    # limit = limit
    # Average hydrogen bond: 1.5-2.5 Angs, not including X-H covalent bond
    # Average electrostatic interaction: < 4 Angs
    # For each system...
    for p in range(len(iad)):
        mol = iad[p]
        # Translate distances into mean contact frequency, with contact defined as when distance < limit
        binary_contacts = distances[mol] < limit
        mean_contacts = np.mean(binary_contacts, 0)
        # Derive coordinates arrays from atom:atom pairs
        x, y = np.transpose(coordinates[mol])[0], np.transpose(coordinates[mol])[1]
        # Plot graph of atom vs. atom, colour coded for (contact frequency)^2
        hm = axis[p].scatter(y, x, c=mean_contacts**2, cmap=cm, vmin=0, vmax=1, marker='s', s=marker_size)
        axis[p].set_title(mol,fontsize=18)
        axis[p].set_xlabel(xlabel,fontsize=14)
        axis[p].set_ylabel(ylabel,fontsize=14)
    cbar = m.colorbar(hm, ax=axis)
    cbar.set_label(colourlabel, fontsize = 14)
    plt.show()
    if saveFig:
        m.savefig(figure_path + mol + '/{sf}_heatmap.png'.format(sf=save_folder))

###############################################################################
#                               ROLAVGPLOT
###############################################################################

def rolAvgPlot2 (data, window_size, no_of_std, ylims, ylabel, legend, figure_name, save): 
    keys = list(data.keys())
    if num_pep == 9:
        g, axs = plt.subplots(3, 3, figsize=(30, 15))
        # reorder = [6,1,5,4,0,2,8,3,7]
        reorder = [0,1,2,3,4,5,6,7,8] 
    elif num_pep == 12:
        g, axs = plt.subplots(3, 4, figsize=(30, 15))
        reorder = [11, 7, 3, 4, 10, 1, 0, 2, 8, 5, 6, 9] 
    axs = axs.ravel()
    # x = np.linspace(0,500000.0,num=25001)
    x = np.linspace(0,1000000.0,num=10001)
    colours=['xkcd:red','xkcd:orange','xkcd:vibrant purple','xkcd:maroon',
            'xkcd:cerulean','xkcd:deep magenta','xkcd:teal','xkcd:green',
             'xkcd:purple','xkcd:grapefruit','xkcd:forest green','xkcd:indigo']
    shade = 0.4    
    c=0
    for i in range(len(keys)):
        nm = keys[i]
        for p in data[nm]:              
            b = reorder[p] 
            df = pd.DataFrame(data[nm][b])
            rolling_mean = df[:].rolling(window_size, center=True).mean()
            rolling_std  = df[:].rolling(window_size, center=True).std()
            num = np.shape(x)[0]
            mean = np.reshape(rolling_mean.values,num)
            minim = np.reshape((rolling_mean + (rolling_std * no_of_std)).values,num)
            maxim = np.reshape((rolling_mean - (rolling_std * no_of_std)).values, num)  
            axs[p].fill_between([a/1000 for a in x],minim, maxim, alpha=(shade/3), facecolor='{}'.format(colours[c]))
            axs[p].plot([a/1000 for a in x],mean, alpha = shade, linewidth=2.0, color='{}'.format(colours[c]), zorder = 20)
            axs[p].set_title("Peptide " + str(b+1), fontsize=22)
            # axs[p].set_xlabel("Simulation Time / $\mu$ s")
            axs[p].set_xlabel("Simulation Time (ns)", fontsize=16)
            axs[p].set_ylabel(ylabel, fontsize=16)
            axs[p].set_ylim(ylims)
            axs[p].set_xlim(0.0, 500)
            axs[p].set_xlim(0.0, 1000)
            axs[p].tick_params(labelsize=12)
            # axs[p].legend(legend, loc=1)  
            c+=1
        c=0
        shade = 1
    plt.subplots_adjust(hspace=0.4)
    if save == True:
        g.savefig(figure_path + mol + '/' + figure_name + ".png")

###############################################################################
#                               PROCESSING HELICITY
###############################################################################

def processHelicity (srcdir, mol,frames):
    print("Processing Helicity")
    col = ['pep'+str(x+1) for x in range(num_pep)]
    df = pd.DataFrame(columns=col)
    print(df.head())
    xtc = '{}{m}/06-MD/{m}_final.xtc'.format(srcdir, m=mol)
    topol = '{}{m}/06-MD/{m}_protein.gro'.format(srcdir, m=mol)

    for chunk in md.iterload(xtc, 100, top=topol):
        in_data = np.empty((100,9))
        for i in range(num_pep): 
            peptide_chunk = chunk.atom_slice(chunk.top.select("resid {} to {}".format(i*pep_length, (i+1)*pep_length-1)))
            dssp = md.compute_dssp(peptide_chunk, simplified=True)
            helicity = (dssp =='H') | (dssp == 'G') | (dssp == 'I')
            in_data[:,i] = np.sum(helicity, axis=1)/1.0 
        df = df.append(pd.DataFrame(in_data, columns=col), ignore_index=True)
    df['tot_hel'] = df.sum(axis=1)
    print(df.head())
    print(df.shape)
    df.iloc[:frames,:].to_csv('./hel_data.csv')
 

def processBonds(limit):
    print("Processing Hydrogen Bonds")
    arcdir = srcdir+'ARC_files/'
    # Directoy names
    sbg = ['EK5_helix']
    # mode = 'salt_bridges'
    mode = 'h_bonds'
    # Defines interval list from which file names are made
    # segments = np.linspace(0,450,10)
    segments = np.linspace(0,950,20)

    distances = {'EK5_helix': np.array([]),}
               # 'E4K4_helix': np.array([]),
                #'mag_helix': np.array([]),
               # }
    coordinates = cp.deepcopy(distances)
    # For each system...
    for r in range(len(sbg)): 
        mol = sbg[r] 
        # For each interval...
        for i in range(len(segments)):  
            lower, upper = int(segments[i])+1, int(segments[i])+50
            target = "{m}_{l}-{u}_{z}.npy".format(m=mol, l=lower, u=upper, z=mode)  
    #        chunk = np.load("{s}{m}/h_bonds/{t}".format(s=srcdir, m=mol, t=target))
            chunk = np.load("{s}{m}/{t}".format(s=arcdir, m=mol, t=target, z=mode))
            # Add distance/time matrix to dictionary entry for that system
            if len(distances[mol]) == 0:
                distances[mol] = chunk
            else:
                distances[mol] = np.append(distances[mol], chunk, axis=0)
            # Set atom:atom pairs (should be the same for every file describing a given system)
            coordinates[mol] = np.load('{s}{m}/{m}_hyb_atom_order.npy'.format(s=arcdir, m=mol, z=mode))
            coordinates[mol] = np.subtract(coordinates[mol], -1)
    binary_contacts = distances[mol] < limit
    df = pd.DataFrame(np.sum(binary_contacts, 1), columns=['tot_hyb'])
    df.to_csv('./hyb_data.csv')

###############################################################################
#                               PLOTTING & RUNNING
###############################################################################
processing = [False,False]
if processing[0]: processHelicity(srcdir,mol,50000)
elif processing[1]: processBonds(0.28)

else:
    helicity_data = pd.read_csv('./hel_data.csv', usecols=['tot_hel'])
    hbond_data = pd.read_csv('./hyb_data.csv', usecols=['tot_hyb'])
    data = helicity_data.merge(hbond_data,left_index=True,right_index=True)
    print(data.head())
    print(data.shape)

    #line_plot = data.plot(y=['tot_hel','tot_hyb'],use_index=True, secondary_y=True)
    #fig = line_plot.get_figure()
    #fig.savefig(figure_path+'test3.png')

    data

'''

#################################################################
helicity = dssp.copy()
#helicity_per_res = dssp.copy()

for mol in helicity:
    for p in helicity[mol]:
        helicity[mol][p] = (dssp[mol][p] =='H')
    #helicity_per_res[mol] = np.sum(helicity[mol], axis=0)[:traj.top.n_residues]/float(traj.n_frames)
        helicity[mol][p] = np.sum(helicity[mol][p], axis=1)/1.0

#rolAvgPlot2(helicity, 200, 1.0, [0,15],"DSSP Score", ["XYZ","Z"], "{}_DSSP".format(mol), True)
#rolAvgPlot2(data, window_size, no_std, ylims, ylabel, figure_title, saveFig)

'''


#plt.plot(EK5_hyb)
#plt.savefig(figure_path+'test2.png')
#############################################################################
'''
system = 'EK5_helix'

EK5_DSSP_overT = []

for i in range(len(helicity[system])):
    EK5_DSSP_overT.append(helicity[system][i])
    
EK5_DSSP_overT = np.array(EK5_DSSP_overT)

EK5_DSSP = np.sum(EK5_DSSP_overT, 0)
'''
