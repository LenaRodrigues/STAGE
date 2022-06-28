import os
import glob
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, SparsePCA, KernelPCA

'''files = glob.glob("/home/lrodrigues/STAGE/DATAclinic/*.clin.merged.picked.txt")
clinicallist = []
for file in files:
    file = pd.read_csv(file, sep="\t")
    file.columns = [list(file)[0]] + [f[:15] for f in list(file)[1:]]
    file = file.T.reset_index()
    file.columns = file.iloc[0, 0:]
    file = file.iloc[1:, :].reset_index(drop=True)
    file.index = file['Hybridization REF']
    file = file.fillna(-476)
    i = 0
    for x in ['vital_status']:
        if x == 1:
            if int(file['days_to_death'][i]) > 365:
                file['vital_status'][i] = 0
        else:
            if int(file['days_to_last_followup'][i]) < 365:
                file = file.drop(file.index[i], axis=0)
    clinicallist.append(file)
df = pd.concat(clinicallist, axis=0)
df.to_csv('/home/lrodrigues/STAGE/DATAclinic/finalclinic.csv', index=False, sep=' ')'''



mRNAseq     = pd.read_csv('/home/lrodrigues/STAGE/DATAmRseq/finalmRNA.csv', low_memory=False)
mRNAseq     = mRNAseq.rename(columns={'gene':'Hybridization REF'})
mRNAseq['Hybridization REF'] = mRNAseq['Hybridization REF'].apply(lambda x: x.lower()[:-3])

RPPA        = pd.read_csv('/home/lrodrigues/STAGE/datarppa/finalrppa.csv', low_memory=False)
RPPA        = RPPA.rename(columns={'Composite.Element.REF':'Hybridization REF'})
RPPA['Hybridization REF'] = RPPA['Hybridization REF'].apply(lambda x: x.lower()[:-3])

methylation = pd.read_csv('/home/lrodrigues/STAGE/DATAMethyl/finalmethyl.csv',low_memory=False)
methylation['Hybridization REF'] = methylation['Hybridization REF'].apply(lambda x: x.lower()[:-3])

miRNAseq    = pd.read_csv('/home/lrodrigues/STAGE/DATAmiRseq/finalmiRNAseq.csv', low_memory=False)
miRNAseq     = miRNAseq.rename(columns={'gene':'Hybridization REF'})
miRNAseq['Hybridization REF'] = miRNAseq['Hybridization REF'].apply(lambda x: x.lower()[:-3])

mRNAseq      = mRNAseq.drop_duplicates(subset=['Hybridization REF'])
RPPA         = RPPA.drop_duplicates(subset=['Hybridization REF'])
methylation  = methylation.drop_duplicates(subset=['Hybridization REF'])
miRNAseq     = miRNAseq.drop_duplicates(subset=['Hybridization REF'])


tmp_list    = np.asarray(list(mRNAseq))
mRNAseq     = mRNAseq[tmp_list[mRNAseq.isna().sum(axis=0) == 0]]

tmp_list = np.asarray(list(RPPA))
RPPA     = RPPA[tmp_list[RPPA.isna().sum(axis=0) == 0]]

tmp_list    = np.asarray(list(methylation))
methylation = methylation[tmp_list[methylation.isna().sum(axis=0) == 0]]

tmp_list    = np.asarray(list(miRNAseq))
miRNAseq    = miRNAseq[tmp_list[miRNAseq.isna().sum(axis=0) == 0]]

label = pd.read_csv('/home/lrodrigues/STAGE/Finalclinic2.csv', header=1, on_bad_lines='skip', low_memory=False)
'''label = label.rename(index ={'HRF': 'Hybridization REF'})
label = label.sort_values(['Hybridization REF', 'Hybridization REF'], axis=0).reset_index(drop=True)
label = label[label['Hybridization REF'].apply(lambda x: 'tcga' in x)].drop_duplicates(subset=['Hybridization REF'], keep ='last').reset_index(drop=True)

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
label.loc[label['days_to_death'].astype(float) <= 5*365, '5yr-mortality'] = 1. '''

'''for view in ['RPPA', 'miRNAseq', 'Methylation', 'mRNAseq']:
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
    df_pca.to_csv('/home/lrodrigues/STAGE/{}_kpca.csv'.format(view), index=False)'''
    
# from sklearn.decomposition import PCA, SparsePCA, KernelPCA

# for view in ['RPPA', 'miRNAseq', 'Methylation', 'mRNAseq']:
#     print(view)
#     if view == 'mRNAseq':
#         df    = mRNAseq.copy(deep=True)
#     elif view == 'miRNAseq':
#         df    = miRNAseq.copy(deep=True)
#     elif view == 'Methylation':
#         df    = methylation.copy(deep=True)
#     elif view == 'RPPA':
#         df    = RPPA.copy(deep=True)

#     z_dim = 100

#     pca   = PCA(n_components=z_dim, random_state=1234)
#     z     =  pca.fit_transform(np.asarray(df.iloc[:, 1:]))

#     df_pca = pd.DataFrame(z, index=df['Hybridization REF']).reset_index()
#     df_pca.to_csv('/home/lrodrigues/STAGE/{}_pca.csv'.format(view), index=False)
    
# from sklearn.decomposition import PCA, SparsePCA, KernelPCA

# for view in ['RPPA', 'miRNAseq', 'Methylation', 'mRNAseq']:
#     print(view)
#     if view == 'mRNAseq':
#         df    = mRNAseq.copy(deep=True)
#     elif view == 'miRNAseq':
#         df    = miRNAseq.copy(deep=True)
#     elif view == 'Methylation':
#         df    = methylation.copy(deep=True)
#     elif view == 'RPPA':
#         df    = RPPA.copy(deep=True)

#     z_dim = 100

#     pca   = SparsePCA(n_components=z_dim, random_state=1234)
#     z     =  pca.fit_transform(np.asarray(df.iloc[:, 1:]))

#     df_pca = pd.DataFrame(z, index=df['Hybridization REF']).reset_index()
#     df_pca.to_csv('/home/lrodrigues/STAGE/{}_spca.csv'.format(view), index=False)

view = 'mRNAseq'
df_pca1  = pd.read_csv('/home/lrodrigues/STAGE/{}_kpca.csv'.format(view), low_memory=False)

view = 'Methylation'
df_pca2  = pd.read_csv('/home/lrodrigues/STAGE/{}_kpca.csv'.format(view), low_memory=False)

view = 'miRNAseq'
df_pca3  = pd.read_csv('/home/lrodrigues/STAGE/{}_kpca.csv'.format(view), low_memory=False)

view = 'RPPA'
df_pca4  = pd.read_csv('/home/lrodrigues/STAGE/{}_kpca.csv'.format(view), low_memory=False)


idx_list_y = label.loc[label['1yr-mortality'] != -1, ['Hybridization REF']]

idx_list1 = df_pca1['Hybridization REF']
idx_list2 = df_pca2['Hybridization REF']
idx_list3 = df_pca3['Hybridization REF']
idx_list4 = df_pca4['Hybridization REF']

idx_list_x = np.unique(idx_list1.tolist() + idx_list2.tolist() + idx_list3.tolist() + idx_list4.tolist())

idx_list     = np.intersect1d(idx_list_x, idx_list_y)
df           = pd.DataFrame(idx_list, columns=['Hybridization REF'])  ##superset of samples that has at least one omics available.


df1 = pd.merge(df, df_pca1, how='left', on='Hybridization REF')
df2 = pd.merge(df, df_pca2, how='left', on='Hybridization REF')
df3 = pd.merge(df, df_pca3, how='left', on='Hybridization REF')
df4 = pd.merge(df, df_pca4, how='left', on='Hybridization REF')
dfy = pd.merge(df, label['Hybridization REF'], how='left', on='Hybridization REF')

np.savez(
    '/home/lrodrigues/STAGE/multi_omics_1yr_mortality.npz',
    mRNAseq     = np.asarray(df1.iloc[:, 1:]),
    Methylation = np.asarray(df2.iloc[:, 1:]),
    miRNAseq    = np.asarray(df3.iloc[:, 1:]),
    RPPA        = np.asarray(df4.iloc[:, 1:]),
    label       = np.asarray(dfy.iloc[:, 1:])
)
