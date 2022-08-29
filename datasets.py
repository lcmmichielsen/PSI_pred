# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 12:05:09 2022

@author: lcmmichielsen
"""

## Script to add the different dataset classes

import pandas as pd
import numpy as np
from torch.utils.data import Dataset
from sklearn.preprocessing import StandardScaler


class EmbOnly(Dataset):
    """Using the SpliceAI embeddings as input to a network only."""
    
    def __init__(self, fn_emb, fn_PSI, idx, cols_PSI=[2,3,4,5,6], 
                 scale=False, train_idx=None):
        """Initialize dataset
        
        Parameters
        ----------
        fn_emb: str
            Name of file embeddings are stored
        fn_PSI: str
            Name of file PSI values are stored
        idx: list
            Indices of dataset used in this fold
        cols_PSI: list, optional
            Indices of columns to use. Default corresponds to PSI values
            for astro, oligo, exc, inh, prog (default = [2,3,4,5,6])
        scale: boolean, optional
            Whether to standard scale the data. If doing so, we also need
            the indices of the training data to fit the scaler (default = False)
        train_idx: list, optional
            If scaling the input data, the training data is needed to train the
            scaler.
        
        """
        
        emb_test = pd.read_csv(fn_emb, index_col=0).values[idx]
        
        if scale:
            if (train_idx == None).all():
                print('Scaling not possible, training idx not defined --> skipping')
            else:
                emb_train = pd.read_csv(fn_emb, index_col=0).values[train_idx]
                scaler = StandardScaler()
                scaler.fit(emb_train)
                emb_test = scaler.transform(emb_test)
        
        self.emb = emb_test
        PSI = pd.read_csv(fn_PSI, index_col=None, usecols=cols_PSI).values[idx]
        tokeep = np.sum(np.isnan(PSI) == False, axis=1) != 0
        self.PSI = PSI[tokeep]
        self.ID = idx     
        
    def __len__(self):
        return len(self.PSI)
        
    def __getitem__(self, index):
        emb = self.emb[index].astype(np.float32)
        PSI = self.PSI[index].astype(np.float32)
        ID = self.ID[index].astype(np.int)
        
        sample = {"emb": emb, "PSI": PSI, "ID": ID}
        return sample