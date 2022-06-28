import numpy as np
import pandas as pd
import random



## all samples MUST include at least one view.
def import_incomplete_handwritten():
    npz = np.load('data_with_missingviews.npz', allow_pickle=True)

    X_set    = npz['X_set'].tolist()
    Y_onehot = npz['Y_onehot']

    M        = len(X_set)

    ### Construct Mask Vector to indicate available (m=1) or missing (m=0) values
    Mask     = np.ones([np.shape(X_set[0])[0], M])
    for m_idx in range(M):
        Mask[np.isnan(X_set[m_idx]).all(axis=1), m_idx] = 0
        X_set[m_idx][Mask[:, m_idx] == 0] = np.mean(X_set[m_idx][Mask[:, m_idx] == 1], axis=0)
    
    return X_set, Y_onehot, Mask



def import_dataset_TCGA(year=1):
    filename = '/home/lrodrigues/STAGE/multi_omics_{}yr_mortality.npz'.format(int(year))
    npz = np.load(filename)
    print(npz.files)
    Y = npz['label']
    Mask = np.zeros((Y.shape[0], 4))
    i = 0
    for m in ['Methylation', 'miRNAseq', 'mRNAseq', 'RPPA']:
        vec = np.any(np.isnan(npz[m]), axis=1)
        Mask[:, i] = vec
        i = i+1

    X_set = {}
    for m in ['Methylation', 'miRNAseq', 'mRNAseq', 'RPPA']:
        tmp = npz['Methylation'.format(1)]
        tmp[np.isnan(tmp[:, 0]), :] = np.nanmean(tmp, axis=0)
        X_set[m] = tmp
    X_set_incomp = {}
    X_set_comp = {}
    for m in ['Methylation', 'miRNAseq', 'mRNAseq', 'RPPA']:
        X_set_comp[m] = X_set[m][np.sum(Mask, axis=1) == 4]
        X_set_incomp[m] = X_set[m][np.sum(Mask, axis=1) != 4]
    Y_comp = Y[np.sum(Mask, axis=1) == 4]
    Y_incomp = Y[np.sum(Mask, axis=1) != 4]
    Mask_comp = Mask[np.sum(Mask, axis=1) == 4]
    Mask_incomp = Mask[np.sum(Mask, axis=1) != 4]
    Y_onehot_incomp = np.zeros([np.shape(Y_incomp)[0], 2])
    Y_onehot_comp = np.zeros([np.shape(Y_comp)[0], 2])

    Y_onehot_incomp[Y_incomp == 0] = 1
    Y_onehot_incomp[Y_incomp == 1] = 1
    Y_onehot_comp[Y_comp == 0] = 1
    Y_onehot_comp[Y_comp == 1] = 1

    return X_set_comp, Y_onehot_comp, Mask_comp, X_set_incomp, Y_onehot_incomp, Mask_incomp
