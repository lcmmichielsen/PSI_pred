# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 12:06:00 2022

@author: lcmmichielsen
"""

import torch

## Script with custom evaluation functions
## Add some metrics that emphasize the differences here?

def loss_difference(out, PSI, current_loss, criterion):

    for i in range(5):
        for j in range(i,5):
            if i != j:
                yt1 = PSI[:,i]
                yt2 = PSI[:,j]
                yp1 = out[:,i]
                yp2 = out[:,j]
    
                notnan = (torch.isnan(yt1) == False) & (torch.isnan(yt2) == False) 
                if torch.sum(notnan):
                    current_loss += criterion(yt1[notnan]-yt2[notnan], yp1[notnan]-yp2[notnan])/5
                    
    return current_loss
