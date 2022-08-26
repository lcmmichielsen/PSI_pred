# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 12:05:09 2022

@author: lcmmichielsen
"""

## Script to add the different dataset classes

import pandas as pd
import numpy as np
from torch.utils.data import Dataset


class emb_only(Dataset):
    """Using the SpliceAI embeddings as input to a network only."""
    
    def __init__(self, fn_emb, fn_PSI, idx, cols=[2,3,4,5,6]):
        """Initialize dataset
        
        Parameters
        ----------
        fn_emb: str
            Name of file embeddings are stored
        fn_PSI: str
            Name of file PSI values are stored
        idx: list
            Indices of dataset used in this fold
        cols: list, optional
            Indices of columns to use. Default corresponds to PSI values
            for astro, oligo, exc, inh, prog (default = [2,3,4,5,6])
        """
        
        self.emb = pd.read_csv(fn_emb, index_col=0).values[idx]
        self.PSI = pd.read_csv(fn_PSI, index_col=None, usecols=cols).values[idx]      
        
    def __len__(self):
        return len(self.PSI)
        
    def __getitem__(self, index):
        emb = self.emb[index].astype(np.float32)
        PSI = self.PSI[index].astype(np.float32)
        
        sample = {"emb": emb, "PSI": PSI}
        return sample