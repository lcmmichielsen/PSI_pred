# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 12:06:26 2022

@author: lcmmichielsen
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
import torch
from torch.utils.data import DataLoader

from datasets import EmbOnly
from networks import MLP
from models import train, test

# Select device
device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
print("[*] Selected device: ", device)

def str2bool(value):
    return value.lower() == 'true'

parser = argparse.ArgumentParser(description='')

parser.add_argument('--fn_emb',         dest='fn_emb',              type=str,       default='/athena/tilgnerlab/store/lim4020/HumanBrainHippProject/downsampled_emb.csv',
                    help='File with the embeddings of spliceAI')
parser.add_argument('--fn_PSI',         dest='fn_PSI',              type=str,       default='/athena/tilgnerlab/store/lim4020/HumanBrainHippProject/downsampled_PSI.csv',
                    help='Directory to save the model and predictions')
parser.add_argument('--cols_PSI',       dest='cols_PSI',            default='2,3,4,5,6',
                    help='Index of columns used to train the model (separated by commas)')
parser.add_argument('--hidden_MLP',     dest='hidden_MLP',          default='0',
                    help='Dimensions of hidden layers MLP model (separated by commas), 0 means no hidden layers')
parser.add_argument('--output_MLP',     dest='output_MLP',          type=int,       default=5,
                    help='Output dimension of MLP model')
parser.add_argument('--dropout',        dest='dropout',             type=float,     default=0.5,
                    help='Dropout rate of MLP')
parser.add_argument('--empdiff',        dest='empdiff',             type=str2bool,  default='True',
                    help='Whether to emphasize difference when calculating the loss')
parser.add_argument('--var_tokeep',     dest='var_tokeep',          type=str,       default="Cons_high,Cons_low,Fake,High,Low,Medium",
                    help='Which type of exons to use during training')
parser.add_argument('--output_dir',     dest='output_dir',          type=str,       default='output',
                    help='Directory to save the model and predictions')


args = parser.parse_args()

cols_PSI = np.asarray(args.cols_PSI.split(','), dtype=int)
output_MLP = args.output_MLP
dropout = args.dropout
empdiff = args.empdiff
var_tokeep = np.asarray(args.var_tokeep.split(','))
output_dir = args.output_dir

if args.hidden_MLP == '0':
    hidden_MLP = []
else:
    hidden_MLP = np.asarray(args.hidden_MLP.split(','), dtype=int)


## Most parameters are now hard-coded --> write complete arg parse later on.
## E.g. epochs, loss function, data type, network, batchsize... 
## Also hyperparameters of MLP are hardcoded (e.g. %dropouts, batchnormalization)
data = 'EmbOnly'
scale = True
runs = 5
network = 'MLP'
fn_emb = args.fn_emb
fn_PSI = args.fn_PSI


try:
    os.mkdir(output_dir)
except:
    print('Output directory already exists')
    
os.chdir(output_dir)


# Split data in train-val-test set (group by genes, stratify by variability)
genes = pd.read_csv(fn_PSI, usecols=[1]).values
varStatus = pd.read_csv(fn_PSI, usecols=[10]).values

fold=0

cv = StratifiedGroupKFold(n_splits=10, random_state=0, shuffle=True)
for train_val_idxs, test_idxs in cv.split(varStatus, varStatus, genes):
    print('Size test set: ')
    print(len(test_idxs))
    
    cv = StratifiedGroupKFold(n_splits=9, random_state=0, shuffle=True)
    ynew = varStatus[train_val_idxs]
    groupsnew = genes[train_val_idxs]
    for train_idxs, val_idxs in cv.split(ynew, ynew, groupsnew):
        train_idxs = train_val_idxs[train_idxs]
        val_idxs = train_val_idxs[val_idxs]
        print('Size validation set: ')
        print(len(val_idxs))
        print('Size training set: ')
        print(len(train_idxs))
        break
    
    # Filter the train-val-test idxs by var status
    train_idxs = train_idxs[np.squeeze(np.isin(varStatus[train_idxs], var_tokeep))]
    val_idxs = val_idxs[np.squeeze(np.isin(varStatus[val_idxs], var_tokeep))]
    test_idxs = test_idxs[np.squeeze(np.isin(varStatus[test_idxs], var_tokeep))]
    
    print('Size test set: ')
    print(len(test_idxs))
    print('Size validation set: ')
    print(len(val_idxs))
    print('Size training set: ')
    print(len(train_idxs))
    
    # Create dataloaders
    if data == 'EmbOnly':
        train_set = EmbOnly(fn_emb, fn_PSI, train_idxs, cols_PSI, scale=True,
                            train_idx=train_idxs)
        val_set = EmbOnly(fn_emb, fn_PSI, val_idxs, cols_PSI, scale=True,
                          train_idx=train_idxs)
        test_set = EmbOnly(fn_emb, fn_PSI, test_idxs, cols_PSI, scale=True,
                           train_idx=train_idxs)
    else:
        print('This data type does not exist --> quitting')
        break
        
    train_loader = DataLoader(train_set, batch_size=64, shuffle=True, drop_last=True)
    train_loader_eval = DataLoader(train_set, batch_size=1, shuffle=False, drop_last=True)
    val_loader = DataLoader(val_set, batch_size=1, shuffle=False, drop_last=True)
    test_loader = DataLoader(test_set, batch_size=1, shuffle=False, drop_last=True)
    
    # Multiple runs per model --> some stochasticity
    for i in range(runs):
        
        logs_dir = 'logs_' + str(fold) + '_' + str(i)
        ckpt_dir = 'ckpt_' + str(fold) + '_' + str(i)
        
        try:
            os.mkdir(ckpt_dir)
            os.mkdir(logs_dir)
        except:
            print('Dir already exists')
        
        if network == 'MLP':
            net = MLP(hidden_dim=hidden_MLP, output_dim=output_MLP, dropout=dropout).to(device)
        else: 
            print('This network is unknown --> quitting')
            break
        
        # Print architecture network & write to file
        # So we don't forget how we trained the model
        print(net)
        print("[*] Number of model parameters:")
        print(sum(p.numel() for p in net.parameters() if p.requires_grad))
        original_stdout = sys.stdout
        with open(logs_dir + '/model.txt', 'w') as tf:
            sys.stdout = tf
            print(net)
            print("[*] Number of model parameters:")
            print(sum(p.numel() for p in net.parameters() if p.requires_grad))
            sys.stdout = original_stdout    
        
        criterion = torch.nn.MSELoss().to(device)
        num_epochs = 50
        init_lr = 0.001
        lr_sched = True
        
        train(device=device, net=net, criterion=criterion,
              learning_rate=init_lr, lr_sched=lr_sched, num_epochs=num_epochs,
              train_loader=train_loader, train_loader_eval=train_loader_eval, 
              valid_loader=val_loader, ckpt_dir=ckpt_dir, logs_dir=logs_dir,
              emp_diff = empdiff)
        
        model_file = ckpt_dir + '/model_best.pth.tar'
        save_file = logs_dir + '/results_testdata_best.pkl'
        test(device=device, net=net, criterion=criterion, model_file=model_file,
                    test_loader=test_loader, save_file=save_file, emp_diff = empdiff)
        save_file = logs_dir + '/results_valdata_best.pkl'
        test(device=device, net=net, criterion=criterion, model_file=model_file,
                    test_loader=val_loader, save_file=save_file, emp_diff = empdiff)
        save_file = logs_dir + '/results_traindata_best.pkl'
        test(device=device, net=net, criterion=criterion, model_file=model_file,
                    test_loader=train_loader_eval, save_file=save_file, emp_diff = empdiff)
        
    fold += 1
        
