"""
###############################################################################
                                BRUSH ANALYSIS
###############################################################################
"""

import mdtraj as md
import numpy as np
import pandas as pd
import os
import glob
import sys
import pathlib
from itertools import repeat

###############################################################################
#                                   INPUTS
###############################################################################

# Source directory for the files on your local system
SRCDIR = '/home/rhys/Desktop/gpfs/rhys/Project/Peptides/carlos_peptides/fibrils/'

# mol = 'brush_3HIS+GLYx3_4x4_0.0A'
# mol = 'brush_3HIS+GLYx3_Ncapped_4x4_0.0A'
FIGDIR = '{}figures/'.format(SRCDIR)
NUM_RES = 16
NUM_PEP = 9
# default_x = np.linspace(0,2500000.0,num=25001)
# plt.style.use('dark_background')

##############################################################################
#                               PROCESSING HELICITY
##############################################################################


def processHelicity(traj, gro, frames):
    print("{:8}{m}".format('INFO:',
                           m='Processing Helicity'))
    # get external constants
    npep = NUM_PEP
    nres = NUM_RES
    # define indexing for the dataframe:
    # |           pep1          |           pep2          |... <- to npep
    # |  res1  |  res2  |  ...  |  res1  |  res2  |  ...  |... <- to nres
    pep = ['pep'+str(x+1) for x in range(npep)]
    multipep = [x for item in pep for x in repeat(item, nres)]
    res = ['res'+str(x+1) for x in range(nres)]
    multires = []
    for i in np.arange(nres):
        multires.extend(res)
        i += 1
    arrays = [multipep, multires]
    tuples = list(zip(*arrays))
    index = pd.MultiIndex.from_tuples(tuples, names=['Peptide', 'Residue'])
    # make dataframe for output using defined indexing
    data = pd.DataFrame(columns=index)
    # convert Path input to string for mdtraj
    xtc = str(traj)
    topol = str(gro)
    # load traj in chunks of timestamps...
    for chunk in md.iterload(xtc, 100, top=topol):
        df = pd.DataFrame(columns=index)
        # ...and one peptide at a time
        for i in range(npep):
            peptide_chunk = chunk.atom_slice(
                                chunk.top.select("resid {} to {}"
                                                 .format(i*nres,
                                                         (i+1)*nres-1)))
            # compute the DSSP values
            dssp = md.compute_dssp(peptide_chunk, simplified=False)
            # add data to new dataframe - store values per timestamp chunk
            for j in np.arange(nres):
                df['pep'+str(i+1), 'res'+str(j+1)] = pd.Series(dssp[:, j]).fillna(0)
        # add data to output dataframe
        data = data.append(df)
    return data
    #df.iloc[:frames, :].to_csv(mol+'_strand_data.csv')
    # helicity = (dssp == 'H') | (dssp == 'G') | (dssp == 'I')
    #helicity = (dssp == 'E') | (dssp == 'B')
    #in_data[:, i] = np.sum(helicity, axis=1)/1.0 
    #new_rows = pd.DataFrame(columns=col) 

def process_Dist():

    pep_length = 293
    diff = 197

    pairs = np.empty((NUM_PEP, 2))
    for n in np.arange(NUM_PEP):
        pairs[n,:] =  [5+n*pep_length, (5+diff)+n*pep_length]
    return pairs


def processHBonds(mol, limit):
    print("Processing Hydrogen Bonds")
    arcdir = SRCDIR+'ARC_files/'
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

def processMetric(SRCDIR, mol, frames, metric):
    print("Processing "+metric)
    col = ['pep'+str(x+1) for x in range(num_pep)]
    df = pd.DataFrame(columns=col)
    print(df.head())
    xtc = '{}{m}/05-MD/MD_final.xtc'.format(SRCDIR, m=mol)
    topol = '{}{m}/05-MD/{m}_protein.gro'.format(SRCDIR, m=mol)

    for chunk in md.iterload(xtc, 100, top=topol):
        in_data = np.empty((100, num_pep))
        for i in range(num_pep):
            peptide_chunk = chunk.atom_slice(chunk.top.select("resid {} to {}"
                                .format(i*pep_length, (i+1)*pep_length-1)))
            #dssp = md.compute_dssp(peptide_chunk, simplified=False)
            if metric == "RGYR":
                in_data[:, i] = md.compute_rg(peptide_chunk)
            elif metric == "RMSD":
                in_data[:, i] = md.rmsd(peptide_chunk, peptide_chunk[0])
            elif metric == "E2E":
                calphas = peptide_chunk.top.select_atom_indices(selection="alpha")
                pair = (calphas[0], calphas[-1])
                mpair = np.matrix(pair)
                in_data[:, i] = np.concatenate(md.compute_distances(peptide_chunk, mpair))

        df = df.append(pd.DataFrame(in_data, columns=col), ignore_index=True)
    df['tot_hel'] = df.mean(axis=1)
    print(df.head())
    print(df.shape)
    df.iloc[:frames, :].to_csv(mol+'_'+metric+'_data.csv')

###############################################################################
#                               PLOTTING & RUNNING
###############################################################################
def run(processing, mol):
    """ run the analysis and/or plot """
    if processing[0]:
        processHelicity(SRCDIR, mol, 50000)
    elif processing[1]:
        processHBonds(mol, 0.28)

    else:
        helicity_data = pd.read_csv(mol+'_hel_data.csv', usecols=['tot_hel'])
        helicity_data['tot_hel'] = (helicity_data['tot_hel'] / helicity_data['tot_hel'].max())*100
        hbond_data = pd.read_csv(mol+'_hyb_data.csv', usecols=['tot_hyb'])
        #print(hbond_data.head())
        data = helicity_data.merge(hbond_data, left_index=True, right_index=True)
        #print(data.shape)

        fontpath = '/home/rhys/.fonts/iosevka/iosevka-term-regular.ttf'
        prop = font_manager.FontProperties(fname=fontpath)
        matplotlib.rcParams['font.family'] = prop.get_name()
        hyb_col = '#ffc000'
        hel_col = '#002C56'

        window_size = 500
        font_sizes = [40, 36]
        x = np.arange(50000)
        for d in ['hel', 'hyb']:
            data[d+'_rm'] = data['tot_'+d].rolling(window_size, center=True).mean()
            data[d+'_std'] = data['tot_'+d].rolling(window_size, center=True).std()
            data[d+'_min'] = data[d+'_rm'] + data[d+'_std']
            data[d+'_max'] = data[d+'_rm'] - data[d+'_std']
        data.to_csv(mol+'_data_out.csv')
        fig, ax1 = plt.subplots(figsize=(20, 16))

        ax1.plot(data['hel_rm'], color=hel_col, linewidth=4.0, zorder=21)
        ax1.fill_between(x, data['hel_min'], data['hel_max'],
                         alpha=.3, facecolor=hel_col, zorder=20)
        ax1.set_ylabel('Percentage Helicity', fontweight='medium',
                       fontsize=font_sizes[0])
        ax1.set_ylim(-2.0, 100.0)
        ax1.set_xlabel('Time (ns)', fontweight='medium', fontsize=font_sizes[0])
        ax1.set_xlim(0.0, 50000.0)
        plt.yticks(fontweight='medium', fontsize=font_sizes[1])
        ax2 = ax1.twinx()

        ax2.plot(data['hyb_rm'], color=hyb_col, linewidth=4.0)
        ax2.fill_between(x, data['hyb_min'], data['hyb_max'], alpha=.3,
                         facecolor=hyb_col)
        ax2.set_ylabel('No. Hydrogen Bonds', fontweight='medium',
                       fontsize=font_sizes[0])
        ax2.set_ylim(-0.1, 5.0)
        plt.yticks(fontweight='medium', fontsize=font_sizes[1])


        ns = np.linspace(0, 1000, 11, dtype='int')
        ts = np.linspace(0, 50000, 11)
        ax1.set_xticks(ticks=ts)
        ax1.set_xticklabels(labels=ns, fontweight='medium', fontsize=font_sizes[1])
    #    ax1.set_title('$(E_4K_4)_2$ Helical Content and Hydrogen Bond Formation', fontweight='medium',fontsize=font_sizes[0], pad=20)
        ax1.set_title('$(EK)_5$ Helical Content and Hydrogen Bond Formation', fontweight='medium',fontsize=font_sizes[0], pad=20)
        #ax1.set_title('AEAK... Helical Content and Hydrogen Bond Formation', fontweight='medium',fontsize=font_sizes[0], pad=20)

        fig.savefig(FIGDIR+mol+'_double_plot.png', bbox_inches='tight', transparent=True, dpi=300)



    fig2, ax2 = plt.subplots(figsize=(20, 16))
    for d in ['capd', 'uncapd']:
        ax2.plot(data[d+'_rm'], color=col[d], linewidth=4.0, zorder=21)
        ax2.fill_between(x, data[d+'_min'], data[d+'_max'],
                        alpha=.3, facecolor=col[d], zorder=20)
    ax2.set_ylabel(ylabels, fontweight='medium',
                    fontsize=font_sizes[0])
    ax2.set_ylim(ylims)
    ax2.set_xlabel('Time (ns)', fontweight='medium', fontsize=font_sizes[0])
    ax2.set_xlim(0.0, 50000.0)
    ns = np.linspace(0, 1000, 11, dtype='int')
    ts = np.linspace(0, 50000, 11)
    ax2.set_xticks(ticks=ts)
    ax2.set_xticklabels(labels=ns, fontweight='medium', fontsize=font_sizes[1])
    ax2.legend(['N Terminal Capped', 'Uncapped'], fontsize=font_sizes[1])
    plt.yticks(fontweight='medium', fontsize=font_sizes[1])
    fig2.savefig(FIGDIR+uncapped_name+"_COMBIND.png", bbox_inches='tight', transparent=True, dpi=300)

def analyse_trajectories():
    hydrophobics = ['I', 'F', 'W']
    orientations = ['PARA', 'INVERTED']
    for hyd in hydrophobics:
        for ornt in orientations:
            traj = SRCDIR + 'EK-{h}H/{o}/EK-{h}H_h2o/05-MD/EK-{h}H_{o}_cut.xtc'.format(h=hyd, o=ornt)
            gro = SRCDIR + 'EK-{h}H/{o}/EK-{h}H_h2o/05-MD/EK-{h}H_h2o_protein.gro'.format(h=hyd, o=ornt)
            df = processHelicity(traj, gro, 200000)
            df.to_hdf(SRCDIR + 'EK-{h}H_{o}_raw.h5'.format(h=hyd, o=ornt), key='raw')


if __name__ == "__main__":

    #analyse_trajectories()
    print(process_Dist())
'''
    hydrophobics = ['I', 'F', 'W']
    orientations = ['PARA', 'INVERTED']
    for hyd in hydrophobics:
        for ornt in orientations:
            with open('Stats.dat', 'a+') as f:
                f.write('\nEK-{}H  {} \n'.format(hyd, ornt))

            df = pd.read_hdf(SRCDIR + 'EK-{h}H_{o}_raw.h5'.format(h=hyd, o=ornt), key='raw')

            for s in pd.unique(df.values.ravel('K')):
                phobic_total = 0
                philic_total = 0
                for i in np.arange(NUM_PEP):
                    for j in np.arange(NUM_RES):
                        count = df['pep'+str(i+1), 'res'+str(j+1)].str.count(s).sum()
                        if j <= 11:
                            phobic_total += count
                        else:
                            philic_total += count
                with open('Stats.dat', 'a+') as f:
                    f.write('   {}  HPHOBIC: {:5.2f}%  HPHILIC {:5.2f}%\n'.format(s, 100*phobic_total/2160000, 100*philic_total/720000))

'''
    #PROCESSING = [False, False]
    #run(PROCESSING, mol)

    #for mol in ['brush_3HIS+GLYx3_4x4_0.0A', 'brush_3HIS+GLYx3_Ncapped_4x4_0.0A']:
    #for mol in ['brush_3HIS+GLYx7_4x4_0.0A',]:
        #pep_length = 13 if "Ncapped" in mol else 12
        #pep_length = 29 if "Ncapped" in mol else 28
        #print("Using PEP_LENGTH:  "+str(pep_length))
        #processHelicity(SRCDIR, mol, 100000)
        #processMetric(SRCDIR, mol, 100000, "E2E")
        #processMetric(SRCDIR, mol, 100000, "RGYR")
        #processMetric(SRCDIR, mol, 100000, "RMSD")

    #plot_brush_csv(mol+'_strand_data.csv', [-0.2, 5.0], '% Stranda (DSSP)')
    #plot_brush_csv(mol+'_hel_data.csv', [-0.2, 7.0], '% Helicity (DSSP)')

    #plot_indcomb(["brush_3HIS+GLYx7_", "4x4_0.0A_hel_data.csv"], 12,
                  #[-2.0, 100.0], '% Helicity (DSSP)', 1000)
    #plot_indcomb(["brush_3HIS+GLYx7_", "4x4_0.0A_strand_data.csv"], 12,
                  #[-2.0, 100.0], '% Strand (DSSP)', 1000)
    #plot_indcomb(["brush_3HIS+GLYx7_", "4x4_0.0A_RMSD_data.csv"], 12,
                  #[0.0, 0.5], 'RMSD from linear / nm', 1000)
    #plot_indcomb(["brush_3HIS+GLYx7_", "4x4_0.0A_RGYR_data.csv"], 12,
                  #[0.0, 2.0], '$R_{gyr}$ / nm', 1000)
    #plot_indcomb(["brush_3HIS+GLYx7_", "4x4_0.0A_E2E_data.csv"], 12,
                  #[0.0, 5.0], 'End-to-end Distance / nm', 1000)
