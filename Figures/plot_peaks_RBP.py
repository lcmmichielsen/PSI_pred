# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 09:06:04 2023

@author: lcmmichielsen
"""

import numpy as np
from tqdm.notebook import tqdm

import seaborn as sns
from matplotlib import pyplot as plt

def plot_peaks_RBP(RBP, exons_info, peaks, PSI_glia, PSI_neur,
                   lowthres=0.5, highthres=0.5, save=False,
                  folder='RBP_HPC/'):
    
    rng = 400
    bins = 50
    
    bins_RBP = np.zeros((len(exons_info),850))

    for k in tqdm(range(len(exons_info))):

        pk_start = peaks[RBP].loc[exons_info[3].values[k]][:,0]
        pk_end = peaks[RBP].loc[exons_info[3].values[k]][:,1]
        exon_start = int(exons_info[1].iloc[k])
        exon_end = int(exons_info[2].iloc[k])

        if exon_end - exon_start < 50:
            continue

        start_bins=[]
        end_bins=[]

        ## Check 1
        start_bin1 = pk_start - exon_start + rng
        start_bin1_keep = (start_bin1 >= 0)  & (start_bin1 < rng)
        start_bins.extend(start_bin1[start_bin1_keep])

        end_bin1 = pk_end - exon_start + rng
        end_bin1_keep = (end_bin1 >= 0) & (end_bin1 < rng)
        end_bins.extend(end_bin1[end_bin1_keep])

        ## Check 2
        start_bin2 = np.round(np.array((bins/(exon_end-exon_start))*(pk_start-exon_start), dtype=float)) + rng
        start_bin2_keep = (start_bin2 >= 400) & (start_bin2 <= 450)
        start_bins.extend(start_bin2[start_bin2_keep])

        end_bin2 = np.round(np.array((bins/(exon_end-exon_start))*(pk_end-exon_start), dtype=float)) + rng
        end_bin2_keep = (end_bin2 >= 400) & (end_bin2 <= 450)
        end_bins.extend(end_bin2[end_bin2_keep])

        ## Check 3
        start_bin3 = pk_start - exon_end + rng + bins
        start_bin3_keep = (start_bin3 > 450) & (start_bin3 < 850)
        start_bins.extend(start_bin3[start_bin3_keep])

        end_bin3 = pk_end - exon_end + rng + bins
        end_bin3_keep = (end_bin3 > 450) & (end_bin3 < 850)
        end_bins.extend(end_bin3[end_bin3_keep])

        start_bins=np.array(start_bins,dtype=int)
        end_bins=np.array(end_bins,dtype=int)

        if len(end_bins) > 0:

            if len(start_bins) == 0:
                bins_RBP[k,:end_bins[0]] = 1
            elif end_bins[0] < start_bins[0]:
                bins_RBP[k,:end_bins[0]] = 1
                j = 1
            else:
                j = 0

        for i in range(len(start_bins)):
            
            if len(end_bins) == 0:
                bins_RBP[k,start_bins[i]:] = 1
            
            elif j >= len(end_bins):
                bins_RBP[k,start_bins[i]:] = 1
            else:
                bins_RBP[k,start_bins[i]:end_bins[j]] = 1
                j += 1

        if exons_info[4].iloc[k] == '-':
            bins_RBP[k] = np.flip(bins_RBP[k])

    var_exons = np.abs(PSI_glia['0'] - PSI_neur['0']) > 0.25
    cons_exons = np.abs(PSI_glia['0'] - PSI_neur['0']) <= 0.25
    
    bins_RBP_var = bins_RBP[var_exons]
    bins_RBP_cons = bins_RBP[cons_exons]

    ### Glia
    PSI_var = PSI_glia[var_exons]
    PSI_cons = PSI_glia[cons_exons]
    
    bins_RBP_sum_lowPSI_var_glia = np.sum(bins_RBP_var[PSI_var['0'] < lowthres], axis=0)/np.sum(PSI_var['0'] < lowthres)
    bins_RBP_sum_higPSI_var_glia = np.sum(bins_RBP_var[PSI_var['0'] >= highthres], axis=0)/np.sum(PSI_var['0'] >= highthres)    
    bins_RBP_sum_lowPSI_cons = np.sum(bins_RBP_cons[PSI_cons['0'] < 0.5], axis=0)/np.sum(PSI_cons['0'] < 0.5)
    bins_RBP_sum_higPSI_cons = np.sum(bins_RBP_cons[PSI_cons['0'] >= 0.5], axis=0)/np.sum(PSI_cons['0'] >= 0.5)

    ### Neurons
    PSI_var = PSI_neur[var_exons]
    PSI_cons = PSI_neur[cons_exons]
    
    bins_RBP_sum_lowPSI_var_neur = np.sum(bins_RBP_var[PSI_var['0'] < lowthres], axis=0)/np.sum(PSI_var['0'] < lowthres)
    bins_RBP_sum_higPSI_var_neur = np.sum(bins_RBP_var[PSI_var['0'] >= highthres], axis=0)/np.sum(PSI_var['0'] >= highthres)    
    
    lim = np.max([bins_RBP_sum_lowPSI_var_glia,
            bins_RBP_sum_higPSI_var_glia,
            bins_RBP_sum_lowPSI_var_neur,
            bins_RBP_sum_higPSI_var_neur,
            bins_RBP_sum_lowPSI_cons,
            bins_RBP_sum_higPSI_cons])
            
    plt.figure(figsize=(3,2))
    plt.axvline(400, c='lightgrey')
    plt.axvline(450, c='lightgrey')
    sns.lineplot(y=bins_RBP_sum_higPSI_var_glia, 
                 x=np.linspace(start=0,stop=849,num=850), color='#D55E00')
    sns.lineplot(y=bins_RBP_sum_lowPSI_var_glia, 
                 x=np.linspace(start=0,stop=849,num=850), color='#029E73')
    plt.ylim([0, lim+0.01])
    plt.title(RBP + ' (Glia)')
    sns.despine()
    
    if save:
        plt.savefig('Evaluate models/Figures/' + folder + RBP + '_glia.pdf', bbox_inches='tight', dpi=1000)
    plt.show()
    
    
    plt.figure(figsize=(3,2))
    plt.axvline(400, c='lightgrey')
    plt.axvline(450, c='lightgrey')
    sns.lineplot(y=bins_RBP_sum_higPSI_var_neur, 
                 x=np.linspace(start=0,stop=849,num=850), color='#D55E00')
    sns.lineplot(y=bins_RBP_sum_lowPSI_var_neur, 
                 x=np.linspace(start=0,stop=849,num=850), color='#029E73')
    plt.ylim([0, lim+0.01])
    plt.title(RBP + ' (Neurons)')
    sns.despine()
    if save:
        plt.savefig('Evaluate models/Figures/' + folder + RBP + '_neur.pdf', bbox_inches='tight', dpi=1000)
    plt.show()
    
    plt.figure(figsize=(3,2))
    plt.axvline(400, c='lightgrey')
    plt.axvline(450, c='lightgrey')
    sns.lineplot(y=bins_RBP_sum_higPSI_cons, 
                 x=np.linspace(start=0,stop=849,num=850), color='#D55E00')
    sns.lineplot(y=bins_RBP_sum_lowPSI_cons, 
                 x=np.linspace(start=0,stop=849,num=850), color='#029E73')
    plt.ylim([0, lim+0.01])
    plt.title(RBP + ' (Cons. exons)')
    sns.despine()
    if save:
        plt.savefig('Evaluate models/Figures/' + folder + RBP + '_cons.pdf', bbox_inches='tight', dpi=1000)
    plt.show()
    
#     return bins_RBP
