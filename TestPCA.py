import os

import numpy as np
import pandas as pd
tumor_list = [
'ACC',
'BLCA',
'BRCA',
'CESC',
'CHOL',
'COAD',
'COADREAD',
'DLBC',
'ESCA',
'FPPP',
'GBM',
'GBMLGG',
'HNSC',
'KICH',
'KIPAN',
'KIRC',
'KIRP',
'LAML',
'LGG',
'LIHC',
'LUAD',
'LUSC',
'MESO',
'OV',
'PAAD',
'PCPG',
'PRAD',
'READ',
'SARC',
'SKCM',
'STAD',
'STES',
'TGCT',
'THCA',
'THYM',
'UCEC',
'UCS',
'UVM']
## 1. FIND SUPERSET OF RPPA FEATURES
feat_list = {}
for tumor in tumor_list:
    filepath = '/home/lrodrigues/STAGE/datarppa/'.format(tumor)
    filename = '{}.rppa.txt'.format(tumor)

    if os.path.exists(filepath + filename):
        tmp = pd.read_csv(filepath + filename, sep='\t')

        tmp.columns = [list(tmp)[0]] + [f[:15] for f in list(tmp)[1:]]
        tmp         = tmp.T.reset_index()
        tmp.columns = tmp.iloc[0, 0:]
        tmp         = tmp.iloc[1:, :].reset_index(drop=True)
        
        feat_list[tumor] = list(tmp)[1:]
        
        if tumor == 'ACC':
            final_feat_list = feat_list[tumor].copy()
            sup_feat_list   = feat_list[tumor].copy()
        else:
            final_feat_list = np.intersect1d(final_feat_list, feat_list[tumor])
            sup_feat_list  += feat_list[tumor]
sup_feat_list = np.unique(sup_feat_list).tolist()
            

for tumor in tumor_list:
    filepath = '/home/lrodrigues/STAGE/datarppa/'.format(tumor)
    filename = '{}.rppa.txt'.format(tumor)
    
    if os.path.exists(filepath + filename):
        tmp = pd.read_csv(filepath + filename, sep='\t')

        tmp.columns = [list(tmp)[0]] + [f[:15] for f in list(tmp)[1:]]
        tmp         = tmp.T.reset_index()
        tmp.columns = tmp.iloc[0, 0:]
        tmp         = tmp.iloc[1:, :].reset_index(drop=True)
        
        tmp_ = pd.DataFrame([], columns=['Composite.Element.REF'] + sup_feat_list)
        tmp_[['Composite.Element.REF'] + feat_list[tumor]] = tmp[['Composite.Element.REF'] + feat_list[tumor]]
        
        if tumor == 'ACC':
#             final_df = tmp[['gene'] + final_feat_list.tolist()]
            final_df = tmp_
        else:
#             final_df = pd.concat([final_df, tmp[['gene'] + final_feat_list.tolist()]], axis=0)
            final_df = pd.concat([final_df, tmp_], axis=0)
    
final_df = final_df.drop_duplicates(subset=['Composite.Element.REF']).reset_index(drop=True)
final_df.to_csv('./FINAL/RPPA.csv', index=False)
RPPA        = pd.read_csv('./FINAL/RPPA.csv')
RPPA        = RPPA.rename(columns={'Composite.Element.REF':'Hybridization REF'})
RPPA['Hybridization REF'] = RPPA['Hybridization REF'].apply(lambda x: x.lower()[:-3])
RPPA         = RPPA.drop_duplicates(subset=['Hybridization REF'])
tmp_list = np.asarray(list(RPPA))
RPPA     = RPPA[tmp_list[RPPA.isna().sum(axis=0) == 0]]
label = pd.read_csv('./FINAL/clinical_label.csv', header=1)
label = label.sort_values(by='Hybridization REF').reset_index(drop=True)
label = label[label['Hybridization REF'].apply(lambda x: 'tcga' in x)].drop_duplicates(subset=['Hybridization REF'], keep ='last').reset_index(drop=True)
'''
    Some of the patients had shifted columns for some reason.
    Manually corrected these errors.
'''

label.loc[label['days_to_last_followup'] == 'endometrial', 'days_to_last_followup'] = label.loc[label['days_to_last_followup'] == 'endometrial', 'days_to_death']
label.loc[label['days_to_last_followup'] == 'endometrial', 'days_to_death'] = label.loc[label['days_to_last_followup'] == 'endometrial', 'vital_status']
label.loc[label['days_to_last_followup'] == 'endometrial', 'vital_status'] = label.loc[label['days_to_last_followup'] == 'endometrial', 'years_to_birth']

label.loc[label['days_to_last_followup'] == 'other  specify', 'days_to_last_followup'] = label.loc[label['days_to_last_followup'] == 'other  specify', 'days_to_death']
label.loc[label['days_to_last_followup'] == 'other  specify', 'days_to_death'] = label.loc[label['days_to_last_followup'] == 'other  specify', 'vital_status']
label.loc[label['days_to_last_followup'] == 'other  specify', 'vital_status'] = label.loc[label['days_to_last_followup'] == 'other  specify', 'years_to_birth']

label['1yr-mortality'] = -1.
label.loc[label['days_to_last_followup'].astype(float) >= 365, '1yr-mortality'] = 0.
label.loc[label['days_to_death'].astype(float) <= 365, '1yr-mortality'] = 1.

label['3yr-mortality'] = -1.
label.loc[label['days_to_last_followup'].astype(float) >= 3*365, '3yr-mortality'] = 0.
label.loc[label['days_to_death'].astype(float) <= 3*365, '3yr-mortality'] = 1.

label['5yr-mortality'] = -1.
label.loc[label['days_to_last_followup'].astype(float) >= 5*365, '5yr-mortality'] = 0.
label.loc[label['days_to_death'].astype(float) <= 5*365, '5yr-mortality'] = 1.
from sklearn.decomposition import PCA, SparsePCA, KernelPCA

for view in ['RPPA', 'miRNAseq', 'Methylation', 'mRNAseq']:
    print(view)
    if view == 'mRNAseq':
        df    = mRNAseq.copy(deep=True)
    elif view == 'miRNAseq':
        df    = miRNAseq.copy(deep=True)
    elif view == 'Methylation':
        df    = methylation.copy(deep=True)
    elif view == 'RPPA':
        df    = RPPA.copy(deep=True)

    z_dim = 100

    pca   = KernelPCA(kernel='poly', n_components=z_dim, random_state=1234)
    z     =  pca.fit_transform(np.asarray(df.iloc[:, 1:]))

    df_pca = pd.DataFrame(z, index=df['Hybridization REF']).reset_index()
    df_pca.to_csv('./FINAL/cleaned/{}_kpca.csv'.format(view), index=False)
    view = 'RPPA'
df_pca4  = pd.read_csv('./FINAL/cleaned/{}_kpca.csv'.format(view))
idx_list_y = label.loc[label['1yr-mortality'] != -1, 'Hybridization REF']
idx_list4 = df_pca4['Hybridization REF']
idx_list     = np.intersect1d(idx_list_x, idx_list_y)
df           = pd.DataFrame(idx_list, columns=['Hybridization REF'])  ##supers
df4 = pd.merge(df, df_pca4, how='left', on='Hybridization REF')
np.savez(
    './FINAL/multi_omics_1yr_mortality.npz',
    RPPA        = np.asarray(df1.iloc[:, 1:]),
    label       = np.asarray(df1.iloc[:, 1:])
)
