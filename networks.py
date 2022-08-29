# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 12:05:24 2022

@author: lcmmichielsen
"""

## Script to add the different network structures to try
import torch.nn as nn

class MLP(nn.Module):
    """Simple MLP to test first."""
    
    def __init__(self, input_dim=320, hidden_dim=[], output_dim=5):
        """Initialize dataset
        
        Parameters
        ----------
        input_dim: int, optional
            Length of input vector (default = 320)
        hidden_dim: list, optional
            List containing  dimensions of hidden layers, e.g. [32,8]. 
            (default = [])
        output_dim: int, optional
            Length of output vector (default = 5)
        """
        
        super(MLP, self).__init__()
        
        fc_layers=[]
        
        # No hidden layers
        if len(hidden_dim) == 0:
            fc_layers.append(nn.Linear(input_dim, output_dim))
        else:
            fc_layers.append(nn.Linear(input_dim, hidden_dim[0]))
            fc_layers.append(nn.Dropout())
            fc_layers.append(nn.ReLU())
            fc_layers.append(nn.BatchNorm1d(num_features=hidden_dim[0]))
            for i in range(len(hidden_dim)-1):
                fc_layers.append(nn.Linear(hidden_dim[i], hidden_dim[i+1]))
                fc_layers.append(nn.Dropout())
                fc_layers.append(nn.ReLU())
                fc_layers.append(nn.BatchNorm1d(num_features=hidden_dim[i+1]))
            fc_layers.append(nn.Linear(hidden_dim[-1], output_dim))
        
        # Ensure that all output is in the range of 0-1 since we predict prob.
        fc_layers.append(nn.Sigmoid())
        
        self.mlp = nn.Sequential(*fc_layers)
        
    def forward(self, x):
        return self.mlp(x)
