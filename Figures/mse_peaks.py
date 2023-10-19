# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 16:28:37 2023

@author: lcmmichielsen
"""

import copy as cp
import numpy as np
from tqdm.notebook import tqdm

from sklearn.metrics import mean_squared_error as mse


def mse_peaks(RBP, exons_info, peaks, PSI_glia, PSI_neur):
    
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
    
    bins_RBP_sum_lowPSI_var = np.sum(bins_RBP_var[PSI_var['0'] < 0.5], axis=0)/np.sum(PSI_var['0'] < 0.5)
    bins_RBP_sum_higPSI_var = np.sum(bins_RBP_var[PSI_var['0'] >= 0.5], axis=0)/np.sum(PSI_var['0'] >= 0.5)    
    bins_RBP_sum_lowPSI_cons = np.sum(bins_RBP_cons[PSI_cons['0'] < 0.5], axis=0)/np.sum(PSI_cons['0'] < 0.5)
    bins_RBP_sum_higPSI_cons = np.sum(bins_RBP_cons[PSI_cons['0'] >= 0.5], axis=0)/np.sum(PSI_cons['0'] >= 0.5)
    
    low_mse = mse(bins_RBP_sum_lowPSI_var, bins_RBP_sum_lowPSI_cons)
    high_mse = mse(bins_RBP_sum_higPSI_var, bins_RBP_sum_higPSI_cons)
    
    exons_low_var = np.sum(np.sum(bins_RBP_var[PSI_var['0'] < 0.5], axis=1) > 0)
    exons_hig_var = np.sum(np.sum(bins_RBP_var[PSI_var['0'] >= 0.5], axis=1) > 0)
    exons_low_cons = np.sum(np.sum(bins_RBP_cons[PSI_cons['0'] < 0.5], axis=1) > 0)
    exons_hig_cons = np.sum(np.sum(bins_RBP_cons[PSI_cons['0'] >= 0.5], axis=1) > 0)

    peaks_low_var = np.sum(bins_RBP_var[PSI_var['0'] < 0.5])
    peaks_hig_var = np.sum(bins_RBP_var[PSI_var['0'] >= 0.5])
    peaks_low_cons = np.sum(bins_RBP_cons[PSI_cons['0'] < 0.5])
    peaks_hig_cons = np.sum(bins_RBP_cons[PSI_cons['0'] >= 0.5])
    
    res_glia = cp.deepcopy(np.array([low_mse, high_mse, exons_low_var, exons_hig_var, 
                    exons_low_cons, exons_hig_cons, peaks_low_var,
                    peaks_hig_var, peaks_low_cons, peaks_hig_cons]))
    ## Neurons
    PSI_var = PSI_neur[var_exons]
    PSI_cons = PSI_neur[cons_exons]
    
    bins_RBP_sum_lowPSI_var = np.sum(bins_RBP_var[PSI_var['0'] < 0.5], axis=0)/np.sum(PSI_var['0'] < 0.5)
    bins_RBP_sum_higPSI_var = np.sum(bins_RBP_var[PSI_var['0'] >= 0.5], axis=0)/np.sum(PSI_var['0'] >= 0.5)    
    bins_RBP_sum_lowPSI_cons = np.sum(bins_RBP_cons[PSI_cons['0'] < 0.5], axis=0)/np.sum(PSI_cons['0'] < 0.5)
    bins_RBP_sum_higPSI_cons = np.sum(bins_RBP_cons[PSI_cons['0'] >= 0.5], axis=0)/np.sum(PSI_cons['0'] >= 0.5)
    
    low_mse = mse(bins_RBP_sum_lowPSI_var, bins_RBP_sum_lowPSI_cons)
    high_mse = mse(bins_RBP_sum_higPSI_var, bins_RBP_sum_higPSI_cons)
    
    exons_low_var = np.sum(np.sum(bins_RBP_var[PSI_var['0'] < 0.5], axis=1) > 0)
    exons_hig_var = np.sum(np.sum(bins_RBP_var[PSI_var['0'] >= 0.5], axis=1) > 0)
    exons_low_cons = np.sum(np.sum(bins_RBP_cons[PSI_cons['0'] < 0.5], axis=1) > 0)
    exons_hig_cons = np.sum(np.sum(bins_RBP_cons[PSI_cons['0'] >= 0.5], axis=1) > 0)

    peaks_low_var = np.sum(bins_RBP_var[PSI_var['0'] < 0.5])
    peaks_hig_var = np.sum(bins_RBP_var[PSI_var['0'] >= 0.5])
    peaks_low_cons = np.sum(bins_RBP_cons[PSI_cons['0'] < 0.5])
    peaks_hig_cons = np.sum(bins_RBP_cons[PSI_cons['0'] >= 0.5])
    
    res_neur = cp.deepcopy(np.array([low_mse, high_mse, exons_low_var, exons_hig_var, 
                    exons_low_cons, exons_hig_cons, peaks_low_var,
                    peaks_hig_var, peaks_low_cons, peaks_hig_cons]))

    
    return res_glia, res_neur