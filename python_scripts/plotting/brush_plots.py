import mdtraj as md
import numpy as np
import matplotlib
import matplotlib.font_manager as font_manager
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
plt.style.use('seaborn-poster')

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
# ##############################################################################
 #                               PROCESSING HELICITY
##############################################################################

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
    df.iloc[:frames,:].to_csv(mol+'_hel_data.csv')


def processBonds(mol,limit):
    print("Processing Hydrogen Bonds")
    arcdir = srcdir+'ARC_files/'
    # mode = 'salt_bridges'
    mode = 'h_bonds'
    # Defines interval list from which file names are made 
    segments = np.linspace(0,950,20)
    # loop over segmented data 
    for i in range(len(segments)):
        lower, upper = int(segments[i])+1, int(segments[i])+50
        target = "{m}_{l}-{u}_{z}.npy".format(m=mol, l=lower, u=upper, z=mode) 
        print('Processing: '+target)
        chunk = np.load("{s}{m}/{t}".format(s=arcdir, m=mol, t=target, z=mode))

        binary_contacts = chunk < limit
        df = pd.DataFrame(np.sum(binary_contacts, 1), columns=['tot_hyb'])
        if i==0: 
            df.to_csv(mol+'_hyb_data.csv', header=True, index=False)
        else:
            df.to_csv(mol+'_hyb_data.csv', mode='a',header=None, index=False)

###############################################################################
#                               PLOTTING & RUNNING
###############################################################################
processing = [False,False]
if processing[0]: processHelicity(srcdir,mol,50000)
elif processing[1]: processBonds(mol,0.28)

else:
    helicity_data = pd.read_csv(mol+'_hel_data.csv', usecols=['tot_hel'])
    helicity_data['tot_hel'] = (helicity_data['tot_hel'] / helicity_data['tot_hel'].max())*100
    hbond_data = pd.read_csv(mol+'_hyb_data.csv', usecols=['tot_hyb'])
    #print(hbond_data.head())
    data = helicity_data.merge(hbond_data,left_index=True,right_index=True)
    #print(data.shape) 

    fontpath = '/home/rhys/.fonts/iosevka/iosevka-term-regular.ttf'
    prop = font_manager.FontProperties(fname=fontpath)
    matplotlib.rcParams['font.family'] = prop.get_name()

    window_size = 500
    font_sizes=[34,30]
    x = np.arange(50000)
    for d in ['hel','hyb']:
        data[d+'_rm'] = data['tot_'+d].rolling(window_size, center=True).mean()
        data[d+'_std'] = data['tot_'+d].rolling(window_size, center=True).std()
        data[d+'_min'] = data[d+'_rm'] + data[d+'_std']
        data[d+'_max'] = data[d+'_rm'] - data[d+'_std']
    data.to_csv(mol+'_data_out.csv')
    fig, ax2 = plt.subplots(figsize=(20,16))

    
    ax2.plot(data['hyb_rm'], color='#ffc000', linewidth=4.0)
    ax2.fill_between(x,data['hyb_min'], data['hyb_max'], alpha=.3, facecolor='#ffc000')
    ax2.set_ylabel('No. Hydrogen Bonds',fontweight='medium', fontsize=font_sizes[0])
    ax2.set_ylim(-0.1, 5.0)
    ax2.set_xlabel('Time (ns)',fontweight='medium', fontsize=font_sizes[0])
    ax2.set_xlim(0.0, 50000.0)
    plt.yticks(fontweight='medium', fontsize=font_sizes[1])  
    
    #ax2.set_xlim(0.0, 50000.0)
    ax1 = ax2.twinx()
    ax1.plot(data['hel_rm'], color='xkcd:royal blue', linewidth=4.0, zorder=21)
    ax1.fill_between(x,data['hel_min'], data['hel_max'], alpha=.3, facecolor='xkcd:royal blue', zorder=20)
    ax1.set_ylabel('Total DSSP Score',fontweight='medium', fontsize=font_sizes[0] )
    ax1.set_ylim(-2.0, 100.0)



    ns = np.linspace(0,1000,11, dtype='int')
    ts = np.linspace(0, 50000, 11)
    ax2.set_xticks(ticks=ts )
    ax2.set_xticklabels(labels=ns,fontweight='medium', fontsize=font_sizes[1])
    plt.yticks(fontweight='medium', fontsize=font_sizes[1]) 

    fig.savefig(figure_path+mol+'_double_plot.png', bbox_inches='tight', transparent=True, dpi=300)

'''
#############################################################
 29 helicity = dssp.copy()
 28 #helicity_per_res = dssp.copy()
 27 
 26 for mol in helicity:
 25     for p in helicity[mol]:
 24         helicity[mol][p] = (dssp[mol][p] =='H')
 23     #helicity_per_res[mol] = np.sum(helicity[mol], axis=0)[:traj.top.n_residues]/float(traj.n_frames)
 22         helicity[mol][p] = np.sum(helicity[mol][p], axis=1)/1.0
 21 
 20 #rolAvgPlot2(helicity, 200, 1.0, [0,15],"DSSP Score", ["XYZ","Z"], "{}_DSSP".format(mol), True)
 19 #rolAvgPlot2(data, window_size, no_std, ylims, ylabel, figure_title, saveFig)
 18 

 16 
 15 
 14 #plt.plot(EK5_hyb)
 13 #plt.savefig(figure_path+'test2.png')
 12 #############################################################################

 10 system = 'EK5_helix'
  9 
  8 EK5_DSSP_overT = []
  7 
  6 for i in range(len(helicity[system])):
  5     EK5_DSSP_overT.append(helicity[system][i])
  4     
  3 EK5_DSSP_overT = np.array(EK5_DSSP_overT)
  2 
  1 EK5_DSSP = np.sum(EK5_DSSP_overT, 0)
251 '''
