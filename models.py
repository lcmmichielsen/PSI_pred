# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 12:06:43 2022

@author: lcmmichielsen
"""

## Script with the train functions

import pickle
import torch
import torch.optim as optim
from tensorboardX import SummaryWriter

def load_checkpoint(net, optimizer=None, scheduler=None, filename='model_last.pth.tar'):
    """Load checkpoint of a trained model."""
    
    start_epoch = 0
    try:
        checkpoint = torch.load(filename)
        start_epoch = checkpoint['epoch']
        net.load_state_dict(checkpoint['state_dict'])
        if optimizer is not None:
            optimizer.load_state_dict(checkpoint['optimizer'])
        if scheduler is not None:
            scheduler.load_state_dict(checkpoint['scheduler'])
        print("\n[*] Loaded checkpoint at epoch %d" % start_epoch)
    except:
        print("[!] No checkpoint found, start epoch 0")

    return start_epoch


def evaluate(device, net, criterion, eval_loader):
    """
    Evaluate performance of the model
    
    Some of the PSI values are NaN. When calculating the loss, we ignore
    those positions.
    
    """
    
    # Eval each sample
    net.eval()
    avg_loss = 0.0
    y_true = []
    y_pred = []
    with torch.no_grad():   # set all 'requires_grad' to False
        for data in eval_loader:
            # Get current batch and transfer to device
            PSI = data['PSI'].to(device, dtype=torch.float)
            tokeep = torch.isnan(PSI) == False
            emb = data['emb'].to(device, dtype=torch.float)
            
            # Forward pass
            out = net(emb)
            
            # Only calculate loss for defined values
            current_loss = criterion(out[tokeep].squeeze(), 
                                     PSI[tokeep].squeeze())
            avg_loss += current_loss.item() / len(eval_loader)
            y_true.append(PSI.cpu().numpy().squeeze())
            y_pred.append(out.cpu().numpy().squeeze())

    return avg_loss, y_true, y_pred


def train(device, net, criterion, learning_rate, lr_sched, num_epochs, 
          train_loader, train_loader_eval, valid_loader, ckpt_dir, logs_dir,
          evaluate_train = True, save_step = 10):
    """Train the model. 
    
    We only calculate the loss over the defined PSI values during training. 
    The model is saved every 'save_step' iterations and the best model so
    far is saved.
    """
    
    best_valid_loss=100
    
    logger = SummaryWriter(logs_dir)
    
    optimizer = optim.Adam(net.parameters(), lr = learning_rate)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', 
                                                     factor=0.1, patience=5)

    start_epoch = load_checkpoint(net, optimizer, scheduler, 
                                 filename = ckpt_dir+'/model_last.pth.tar')
    
    # Evaluate validation set before start training
    print("[*] Evaluating epoch %d..." % start_epoch)
    avg_valid_loss, _, _ = evaluate(device, net, criterion, valid_loader)
    print("--- Average valid loss:                  %.4f" % avg_valid_loss)

    # Training epochs
    for epoch in range(start_epoch, num_epochs):
        net.train()
        # Print current learning rate
        print("[*] Epoch %d..." % (epoch + 1))
        for param_group in optimizer.param_groups:
            print('--- Current learning rate: ', param_group['lr'])

        for data in train_loader:
            # Get current batch and transfer to device
            PSI = data['PSI'].to(device, dtype=torch.float)
            tokeep = torch.isnan(PSI) == False
            emb = data['emb'].to(device, dtype=torch.float)

            with torch.set_grad_enabled(True):  # no need to specify 'requires_grad' in tensors
                # Set the parameter gradients to zero
                optimizer.zero_grad()

                # Forward pass
                out = net(emb)
                current_loss = criterion(out[tokeep].squeeze(), 
                                         PSI[tokeep].squeeze())
                # print(current_loss)
                # Backward pass and optimize
                current_loss.backward()
                optimizer.step()

        # Save last model
        state = {'epoch': epoch + 1, 'state_dict': net.state_dict(),
                 'optimizer': optimizer.state_dict(), 'scheduler': scheduler.state_dict()}
        torch.save(state, ckpt_dir + '/model_last.pth.tar')

        # Save model at epoch
        if (epoch + 1) % save_step == 0:
            print("[*] Saving model epoch %d..." % (epoch + 1))
            torch.save(state, ckpt_dir + '/model_epoch%d.pth.tar' % (epoch + 1))
        
        # Evaluate all training set and validation set at epoch
        print("[*] Evaluating epoch %d..." % (epoch + 1))
        if evaluate_train:
            avg_train_loss, _, _ = evaluate(device, net, criterion, train_loader_eval)
            print("--- Average train loss:                  %.4f" % avg_train_loss)
            
            logger.add_scalar('train_loss_epoch', avg_train_loss, epoch + 1)

        avg_valid_loss, _, _ = evaluate(device, net, criterion, valid_loader)
        print("--- Average valid loss:                  %.4f" % avg_valid_loss)
        
        # Check if best model
        if avg_valid_loss < best_valid_loss:
            best_valid_loss = avg_valid_loss
            state = {'epoch': epoch + 1, 'state_dict': net.state_dict(),
                     'optimizer': optimizer.state_dict(), 'scheduler': scheduler.state_dict()}
            torch.save(state, ckpt_dir + '/model_best.pth.tar')

        
        logger.add_scalar('valid_loss_epoch', avg_valid_loss, epoch + 1)
        
        # LR scheduler on plateau (based on validation loss)
        if lr_sched:
            scheduler.step(avg_valid_loss)

    print("[*] Finish training.")
    
    

def test(device, net, criterion, model_file, test_loader, save_file=None):
    """Test performance of the current model and write results to a file."""
    
    # Load pretrained model
    _ = load_checkpoint(net, filename=model_file)
    
    # Evaluate model
    avg_test_loss, y_true, y_pred = evaluate(device, net, criterion, test_loader)

    # Save predictions
    if save_file is not None:
        pickle.dump({'y_true': y_true, 'y_pred': y_pred}, open(save_file, 'wb'))

    # Display evaluation metrics
    print("--- Average test loss:                  %.4f" % avg_test_loss)