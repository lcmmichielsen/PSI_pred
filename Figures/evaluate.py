import os
import copy as cp
import h5py
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import r2_score, mean_squared_error, classification_report, confusion_matrix, f1_score
from scipy.spatial import ConvexHull

pd.set_option('mode.chained_assignment', None)


def evaluate(prefixs, folders, multiheads, heads, files_glia, files_neur,
             runs=5, folds = [0,1,2,3,4,5,6,7,8,9], split='test', 
             RBP_only=False, exon_RBP=[]):
    
    res_all_neurMSE = pd.DataFrame(data=np.zeros((len(folders),10)),
                               index=folders)
    res_all_gliaMSE = pd.DataFrame(data=np.zeros((len(folders),10)),
                               index=folders)
    res_all_neurSP = pd.DataFrame(data=np.zeros((len(folders),10)),
                               index=folders)
    res_all_gliaSP = pd.DataFrame(data=np.zeros((len(folders),10)),
                               index=folders)
    res_var_neurMSE = pd.DataFrame(data=np.zeros((len(folders),10)),
                               index=folders)
    res_var_gliaMSE = pd.DataFrame(data=np.zeros((len(folders),10)),
                               index=folders)
    res_var_neurSP = pd.DataFrame(data=np.zeros((len(folders),10)),
                               index=folders)
    res_var_gliaSP = pd.DataFrame(data=np.zeros((len(folders),10)),
                               index=folders)
    F1_scores = pd.DataFrame(data=np.zeros((len(folders),2)),
                            index=folders)
    
    y_true_glia_all = []
    y_pred_glia_all = []
    
    y_true_neur_all = []
    y_pred_neur_all = []
    
    var_genes = []
    
    
    for i,f in enumerate(folders):
        print(f)
        res_neur = summarize_res(prefix=prefixs[i],
                                folder=f+'/neur/',
                                multihead=multiheads[i], 
                                head=heads[i],
                                runs=runs,
                                folds=folds,
                                split=split,
                                RBP_only=RBP_only,
                                exon_RBP=exon_RBP,
                                PSI_glia=files_glia[i],
                                PSI_neur=files_neur[i])
        
        res_glia = summarize_res(prefix=prefixs[i],
                                folder=f+'/glia/',
                                multihead=multiheads[i], 
                                head=heads[i],
                                runs=runs,
                                folds=folds,
                                split=split,
                                RBP_only=RBP_only,
                                exon_RBP=exon_RBP,
                                PSI_glia=files_glia[i],
                                PSI_neur=files_neur[i])
        
        res_var_neurMSE.iloc[i] = res_neur[0]['MSE'].values
        res_var_gliaMSE.iloc[i] = res_glia[0]['MSE'].values
        res_var_neurSP.iloc[i] = res_neur[0]['Spearman'].values
        res_var_gliaSP.iloc[i] = res_glia[0]['Spearman'].values

        res_all_neurMSE.loc[f] = res_neur[1]['MSE'].values
        res_all_gliaMSE.loc[f] = res_glia[1]['MSE'].values
        res_all_neurSP.loc[f] = res_neur[1]['Spearman'].values
        res_all_gliaSP.loc[f] = res_glia[1]['Spearman'].values
    
        delta_true = (res_neur[2]-res_glia[2])>0
        delta_pred = (res_neur[3]-res_glia[3])>0
        
        F1_scores.iloc[i][0] = f1_score(delta_true[res_neur[5]], delta_pred[res_neur[5]], pos_label=0)
        F1_scores.iloc[i][1] = f1_score(delta_true[res_neur[5]], delta_pred[res_neur[5]], pos_label=1)
        
        y_true_glia_all.append(res_glia[2])
        y_pred_glia_all.append(res_glia[3])
        
        y_true_neur_all.append(res_neur[2])
        y_pred_neur_all.append(res_neur[3])
        
        var_genes.append(res_neur[5])

    try:
        
        y_true_glia_all = pd.DataFrame(np.transpose(y_true_glia_all),
                                       index=res_glia[4], columns=folders)
        y_pred_glia_all = pd.DataFrame(np.transpose(y_pred_glia_all),
                                       index=res_glia[4], columns=folders)
        y_true_neur_all = pd.DataFrame(np.transpose(y_true_neur_all),
                                       index=res_neur[4], columns=folders)
        y_pred_neur_all = pd.DataFrame(np.transpose(y_pred_neur_all),
                                       index=res_neur[4], columns=folders)
        var_genes = pd.DataFrame(np.transpose(var_genes), 
                                 index=res_neur[4], columns=folders)
    except:
        1+1
    
    return res_all_neurMSE, res_all_gliaMSE, res_all_neurSP, res_all_gliaSP, res_var_neurMSE, res_var_gliaMSE, res_var_neurSP, res_var_gliaSP, F1_scores, y_true_neur_all, y_pred_neur_all, y_true_glia_all, y_pred_glia_all, var_genes


def summarize_res(prefix, folder, multihead=False, head='/head0', runs=5,
                 folds = [0,1,2,3,4,5,6,7,8,9], split='test', RBP_only=False, exon_RBP=[], PSI_glia='/PSI_glia_norm.csv',
                 PSI_neur='/PSI_neur_norm.csv'):
    
    glia_true = pd.read_csv(PSI_glia, index_col=0)
    neur_true = pd.read_csv(PSI_neur, index_col=0)
    exon_var = glia_true.index[np.squeeze((np.abs(glia_true-neur_true) > 0.25).values)]
    
    if RBP_only:
        exon_var = np.intersect1d(exon_RBP, exon_var)
    
    true_all = []
    pred_all = []
    genes_all = []
    var_all = []

    metrics_all = pd.DataFrame(data=np.zeros((10,6)),
                              columns=['Pearson','Pearson_best','Spearman','MSE','R2','folder'])
    metrics_var = pd.DataFrame(data=np.zeros((10,6)),
                              columns=['Pearson','Pearson_best','Spearman','MSE','R2','folder'])
    metrics_all = metrics_all.astype({'folder': '<U100'})
    metrics_var = metrics_var.astype({'folder': '<U100'})
    metrics_all.iloc[:]['folder'] = folder
    metrics_var.iloc[:]['folder'] = folder
        
    for fold in folds:
        pred_average = 0
        
        if multihead:
            genes = pd.read_csv(prefix+'/'+folder+'/fold'+str(fold)+'/run'+str(0)+head+'/genes.csv', index_col=0)
        else:
            genes = pd.read_csv(prefix+'/'+folder+'/fold'+str(fold)+'/run'+str(0)+'/genes.csv', index_col=0)
        genes = genes[genes['split'] == split]
        
        try:
            index = pd.DataFrame(genes.index)[0].str.split(pat="\.|_", expand=True)
            genes.index = index[0] + '_' + index[1] + '_' + index[2] + '_' + index[3] + '_' + index[5] 
        except:
            1+1

        idx_var = np.isin(genes.index, exon_var)
        # print(idx_var)
        # print(genes.index)
        # print(exon_var)
        
        if RBP_only:
            idx_RBP = np.isin(genes.index, exon_RBP)

        res_ind_runs = np.zeros((5,1))
        res_ind_runs_var = np.zeros((5,1))

        for run in range(runs):
            
            if multihead:
                f1 = h5py.File(prefix+'/'+folder+'/fold'+str(fold)+'/run'+str(run)+head+'/targets.h5', 'r+')
                f2 = h5py.File(prefix+'/'+folder+'/fold'+str(fold)+'/run'+str(run)+head+'/preds.h5', 'r+')
            else: 
                f1 = h5py.File(prefix+'/'+folder+'/fold'+str(fold)+'/run'+str(run)+'/targets.h5', 'r+')
                f2 = h5py.File(prefix+'/'+folder+'/fold'+str(fold)+'/run'+str(run)+'/preds.h5', 'r+')

            target_genes = f1['genes'][()]
            target_genes = np.array(target_genes, dtype='<U50')
            target_values = np.squeeze(f1['targets'][()])
            f1.close()

            pred_genes = f2['genes'][()]
            pred_genes = np.array(pred_genes, dtype='<U50')
            pred_values = np.squeeze(f2['preds'][()])
            f2.close()
            
#             pred_values = np.clip(pred_values, a_min=0, a_max=1)
            pred_average += pred_values/runs
            res_ind_runs[run], _ = pearsonr(target_values,pred_values)
            res_ind_runs_var[run], _ = pearsonr(target_values[idx_var], pred_values[idx_var])
#             count += 1

#             print(np.var(pred_values))

        # print(pred_genes)

        # assert(np.all(pred_genes == genes.index))
        assert(np.all(pred_genes == target_genes))
        
        if RBP_only:
            
            metrics_all['Pearson'].iloc[fold], _ = pearsonr(target_values[idx_RBP], pred_average[idx_RBP])
            metrics_all['Spearman'].iloc[fold], _ = spearmanr(target_values[idx_RBP], pred_average[idx_RBP])
            metrics_all['R2'].iloc[fold] = r2_score(target_values[idx_RBP], pred_average[idx_RBP])
            metrics_all['MSE'].iloc[fold] = mean_squared_error(target_values[idx_RBP], pred_average[idx_RBP])
            
            true_all.extend(target_values[idx_RBP])
            pred_all.extend(pred_values[idx_RBP])
            genes_all.extend(target_genes[idx_RBP])
            
        else:
            metrics_all['Pearson'].iloc[fold], _ = pearsonr(target_values, pred_average)
            metrics_all['Spearman'].iloc[fold], _ = spearmanr(target_values, pred_average)
            metrics_all['R2'].iloc[fold] = r2_score(target_values, pred_average)
            metrics_all['MSE'].iloc[fold] = mean_squared_error(target_values, pred_average)
            
            true_all.extend(target_values)
            pred_all.extend(pred_values)
            genes_all.extend(target_genes)
        
        metrics_all['Pearson_best'].iloc[fold] = np.max(res_ind_runs)
        
        idx_var = np.isin(genes.index, exon_var)
#         print(np.sum(idx_var))
        metrics_var['Pearson'].iloc[fold], _ = pearsonr(target_values[idx_var], pred_average[idx_var])
        metrics_var['Spearman'].iloc[fold], _ = spearmanr(target_values[idx_var], pred_average[idx_var])
        metrics_var['R2'].iloc[fold] = r2_score(target_values[idx_var], pred_average[idx_var])
        metrics_var['MSE'].iloc[fold] = mean_squared_error(target_values[idx_var], pred_average[idx_var])
        metrics_var['Pearson_best'].iloc[fold] = np.max(res_ind_runs_var)

        if RBP_only:
            
            var_all.extend(idx_var[idx_RBP])
        else:
            var_all.extend(idx_var)
        
    return metrics_var, metrics_all, np.asarray(true_all), np.asarray(pred_all), np.asarray(genes_all), np.asarray(var_all), res_ind_runs_var

def scatterplot(x, y, idx_var=False, title='title', s=5, hue=None,
               xlabel='True PSI', ylabel='Pred PSI', kdeplot=True):
    x = np.asarray(x)
    y = np.asarray(y)
    
    if np.any(idx_var):
    
        plt.figure(figsize=(5,5))
        if np.any(hue != None):
            hue = hue[idx_var]
        sns.scatterplot(x=x[idx_var], y=y[idx_var], s=s, hue=hue)
        if kdeplot:
                sns.kdeplot(x=x[idx_var], y=y[idx_var], cut=0.1, alpha=0.5, fill=True)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.show()
        
    else:
        
        plt.figure(figsize=(5,5))
        sns.scatterplot(x=x, y=y, s=s, hue=hue)
        if kdeplot:
            sns.kdeplot(x=x, y=y, cut=0.1, alpha=0.5, fill=True)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.show()

def reorder_df(df, names):
#     names = ['DL - seq', 
#              'DL - seq+RBP', 
#              'DL - RBP',
#              'LR - all exons',
#              'LR - var. exons > 0.1']
    df = df.rename(columns={0: 'v0', 1: 'v1',
                            2: 'v2', 3: 'v3',
                            4: 'v4', 5: 'v5',
                            6: 'v6', 7: 'v7',
                            8: 'v8', 9: 'v9',})
    df['id'] = df.index
    df['name'] = names
    df = pd.wide_to_long(df, stubnames='v', i='id', j='fold')
    df = df.reset_index()
    
    return df


def plot_results(res_all_neur, res_all_glia, res_var_neur, res_var_glia,
                 res_all_neurR2, res_all_gliaR2, res_var_neurR2, res_var_gliaR2,
                 res_all_neurMSE, res_all_gliaMSE, res_var_neurMSE, res_var_gliaMSE,
                 res_all_neurSP, res_all_gliaSP, res_var_neurSP, res_var_gliaSP, names, order):
    
    res_all_neur2 = cp.deepcopy(res_all_neur)
    res_all_glia2 = cp.deepcopy(res_all_glia)
    res_var_neur2 = cp.deepcopy(res_var_neur)
    res_var_glia2 = cp.deepcopy(res_var_glia)

    res_all_neurR22 = cp.deepcopy(res_all_neurR2)
    res_all_gliaR22 = cp.deepcopy(res_all_gliaR2)
    res_var_neurR22 = cp.deepcopy(res_var_neurR2)
    res_var_gliaR22 = cp.deepcopy(res_var_gliaR2)

    res_all_neurMSE2 = cp.deepcopy(res_all_neurMSE)
    res_all_gliaMSE2 = cp.deepcopy(res_all_gliaMSE)
    res_var_neurMSE2 = cp.deepcopy(res_var_neurMSE)
    res_var_gliaMSE2 = cp.deepcopy(res_var_gliaMSE)

    res_all_neurSP2 = cp.deepcopy(res_all_neurSP)
    res_all_gliaSP2 = cp.deepcopy(res_all_gliaSP)
    res_var_neurSP2 = cp.deepcopy(res_var_neurSP)
    res_var_gliaSP2 = cp.deepcopy(res_var_gliaSP)

    res_all_neur = reorder_df(res_all_neur, names)
    res_all_glia = reorder_df(res_all_glia, names)
    res_var_neur = reorder_df(res_var_neur, names)
    res_var_glia = reorder_df(res_var_glia, names)

    res_all_neurR2 = reorder_df(res_all_neurR2, names)
    res_all_gliaR2 = reorder_df(res_all_gliaR2, names)
    res_var_neurR2 = reorder_df(res_var_neurR2, names)
    res_var_gliaR2 = reorder_df(res_var_gliaR2, names)

    res_all_neurMSE = reorder_df(res_all_neurMSE, names)
    res_all_gliaMSE = reorder_df(res_all_gliaMSE, names)
    res_var_neurMSE = reorder_df(res_var_neurMSE, names)
    res_var_gliaMSE = reorder_df(res_var_gliaMSE, names)

    res_all_neurSP = reorder_df(res_all_neurSP, names)
    res_all_gliaSP = reorder_df(res_all_gliaSP, names)
    res_var_neurSP = reorder_df(res_var_neurSP, names)
    res_var_gliaSP = reorder_df(res_var_gliaSP, names)

    res_all_neur['Cell type'] = 'Neurons'
    res_all_glia['Cell type'] = 'Glia'
    res_var_neur['Cell type'] = 'Neurons'
    res_var_glia['Cell type'] = 'Glia'

    res_all_neurR2['Cell type'] = 'Neurons'
    res_all_gliaR2['Cell type'] = 'Glia'
    res_var_neurR2['Cell type'] = 'Neurons'
    res_var_gliaR2['Cell type'] = 'Glia'

    res_all_neurMSE['Cell type'] = 'Neurons'
    res_all_gliaMSE['Cell type'] = 'Glia'
    res_var_neurMSE['Cell type'] = 'Neurons'
    res_var_gliaMSE['Cell type'] = 'Glia'

    res_all_neurSP['Cell type'] = 'Neurons'
    res_all_gliaSP['Cell type'] = 'Glia'
    res_var_neurSP['Cell type'] = 'Neurons'
    res_var_gliaSP['Cell type'] = 'Glia'

    res_all = pd.concat((res_all_neur, res_all_glia))
    res_var = pd.concat((res_var_neur, res_var_glia))

    res_allR2 = pd.concat((res_all_neurR2, res_all_gliaR2))
    res_varR2 = pd.concat((res_var_neurR2, res_var_gliaR2))

    res_allMSE = pd.concat((res_all_neurMSE, res_all_gliaMSE))
    res_varMSE = pd.concat((res_var_neurMSE, res_var_gliaMSE))

    res_allSP = pd.concat((res_all_neurSP, res_all_gliaSP))
    res_varSP = pd.concat((res_var_neurSP, res_var_gliaSP))

    sns.boxplot(data=res_all, x='v', y='name', hue='Cell type', 
                order = order,
               palette=['#9DC3E6', '#F4B68C'])
    sns.despine()
    plt.title('All exons')
    plt.xlabel('Pearson correlation')
    plt.ylabel('')
    plt.legend(bbox_to_anchor=(0.025, 0.35), loc='upper left', borderaxespad=0)
    plt.show()

    sns.boxplot(data=res_var, x='v', y='name', hue='Cell type',
                order = order,
               palette=['#9DC3E6', '#F4B68C'])
    sns.despine()
    plt.title('Variable exons')
    plt.xlabel('Pearson correlation')
    plt.ylabel('')
    plt.legend([], [], frameon=False)
    plt.show()

    sns.boxplot(data=res_allMSE, x='v', y='name', hue='Cell type', 
                order = order,
               palette=['#9DC3E6', '#F4B68C'])
    sns.despine()
    plt.title('All exons')
    plt.xlabel('MSE')
    plt.ylabel('')
    plt.legend([], [], frameon=False)
    plt.show()

    sns.boxplot(data=res_varMSE, x='v', y='name', hue='Cell type',
                order = order,
               palette=['#9DC3E6', '#F4B68C'])
    sns.despine()
    plt.title('Variable exons')
    plt.xlabel('MSE')
    plt.ylabel('')
    plt.legend([], [], frameon=False)
    plt.show()

    sns.boxplot(data=res_allSP, x='v', y='name', hue='Cell type', 
                order = order,
               palette=['#9DC3E6', '#F4B68C'])
    sns.despine()
    plt.title('All exons')
    plt.xlabel('Spearman correlation')
    plt.ylabel('')
    plt.legend([], [], frameon=False)
    plt.show()

    sns.boxplot(data=res_varSP, x='v', y='name', hue='Cell type',
                order = order,
               palette=['#9DC3E6', '#F4B68C'])
    sns.despine()
    plt.title('Variable exons')
    plt.xlabel('Spearman correlation')
    plt.ylabel('')
    plt.legend([], [], frameon=False)
    plt.show()

    sns.boxplot(data=res_allR2, x='v', y='name', hue='Cell type', 
                order = order,
               palette=['#9DC3E6', '#F4B68C'])
    sns.despine()
    plt.title('All exons')
    plt.xlabel('R2')
    plt.ylabel('')
    plt.legend([], [], frameon=False)
    plt.show()

    sns.boxplot(data=res_varR2, x='v', y='name', hue='Cell type',
                order = order,
               palette=['#9DC3E6', '#F4B68C'])
    sns.despine()
    plt.title('Variable exons')
    plt.xlabel('R2')
    plt.ylabel('')
    plt.legend([], [], frameon=False)
    plt.show()