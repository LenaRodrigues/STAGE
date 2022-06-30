
from class_DeepIMV_AISTATS import DeepIMV_AISTATS
from helper import f_get_minibatch_set, evaluate
import import_data as impt
from sklearn.model_selection import train_test_split
import os
import sys
import random
import tensorflow as tf
import numpy as np
import warnings
warnings.filterwarnings('ignore')

year = 1
DATASET_PATH = 'TCGA_{}YR'.format(int(year))
DATASET = 'TCGA'

X_set_comp, Y_onehot_comp, Mask_comp, X_set_incomp, Y_onehot_incomp, Mask_incomp = impt.import_dataset_TCGA(year)

MODE = 'incomplete'
model_name = 'DeepIMV_AISTATS'

M = len(X_set_comp)

SEED = 1234
OUTITERATION = 5

RESULTS_AUROC_RAND = np.zeros([4, OUTITERATION+2])
RESULTS_AUPRC_RAND = np.zeros([4, OUTITERATION+2])

out_itr = 1

tr_X_set, te_X_set, va_X_set = {}, {}, {}
for m in ['Methylation', 'miRNAseq', 'mRNAseq', 'RPPA']:
    tr_X_set[m], te_X_set[m] = train_test_split(
        X_set_comp[m], test_size=0.2, random_state=SEED + out_itr)
    tr_X_set[m], va_X_set[m] = train_test_split(
        tr_X_set[m], test_size=0.2, random_state=SEED + out_itr)

tr_Y_onehot, te_Y_onehot, tr_M, te_M = train_test_split(
    Y_onehot_comp, Mask_comp, test_size=0.2, random_state=SEED + out_itr)

tr_Y_onehot, va_Y_onehot, tr_M, va_M = train_test_split(
    tr_Y_onehot, tr_M, test_size=0.2, random_state=SEED + out_itr)

if MODE == 'incomplete':
    for m in ['Methylation', 'miRNAseq', 'mRNAseq', 'RPPA']:
        tr_X_set[m] = np.concatenate([tr_X_set[m], X_set_incomp[m]], axis=0)

    tr_Y_onehot = np.concatenate([tr_Y_onehot, Y_onehot_incomp], axis=0)
    tr_M = np.concatenate([tr_M, Mask_incomp], axis=0)

    print(tr_M.shape)
elif MODE == 'complete':
    print(tr_M.shape)
else:
    raise ValueError('WRONG MODE!!!')


save_path = '{}/M{}_{}/{}/'.format(DATASET_PATH, M, MODE, model_name)


if not os.path.exists(save_path + 'itr{}/'.format(out_itr)):
    os.makedirs(save_path + 'itr{}/'.format(out_itr))

(5850, 4)
# training coefficients
alpha = 1.0
beta = 0.01  # IB coefficient
lr_rate = 1e-4
k_prob = 0.7


# network parameters
mb_size = 32
steps_per_batch = int(np.shape(tr_M)[0]/mb_size)
steps_per_batch = 500

x_dim_set = [tr_X_set[m].shape[1] for m in ['Methylation', 'miRNAseq', 'mRNAseq', 'RPPA']]
y_dim = np.shape(tr_Y_onehot)[1]
y_type = 'binary'
z_dim = 100

h_dim_p = 100
num_layers_p = 2

h_dim_e = 300
num_layers_e = 3

input_dims = {
    'x_dim_set': x_dim_set,
    'y_dim': y_dim,
    'y_type': y_type,
    'z_dim': z_dim,

    'steps_per_batch': steps_per_batch
}

network_settings = {
    'h_dim_p1': h_dim_p,
    'num_layers_p1': num_layers_p,  # view-specific
    'h_dim_p2': h_dim_p,
    'num_layers_p2': num_layers_p,  # multi-view
    'h_dim_e': h_dim_e,
    'num_layers_e': num_layers_e,
    'fc_activate_fn': tf.nn.relu,
    'reg_scale': 0.,  # 1e-4,
}

tf.reset_default_graph()

# gpu_options = tf.GPUOptions()
gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=0.22)
sess = tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))

model = DeepIMV_AISTATS(sess, "DeepIMV_AISTATS", input_dims, network_settings)

saver = tf.train.Saver()
sess.run(tf.global_variables_initializer())

saver = tf.train.Saver()
sess.run(tf.global_variables_initializer())

ITERATION = 500000
STEPSIZE = 500

min_loss = 1e+8
max_acc = 0.0
max_flag = 20

tr_avg_Lt, tr_avg_Lp, tr_avg_Lkl, tr_avg_Lps, tr_avg_Lkls, tr_avg_Lc = 0, 0, 0, 0, 0, 0
va_avg_Lt, va_avg_Lp, va_avg_Lkl, va_avg_Lps, va_avg_Lkls, va_avg_Lc = 0, 0, 0, 0, 0, 0

stop_flag = 0
for itr in range(ITERATION):
    x_mb_set, y_mb, m_mb = f_get_minibatch_set(
        mb_size, tr_X_set, tr_Y_onehot, tr_M)
    
    x_mb_set[0]=x_mb_set.pop('Methylation')
    x_mb_set[1]=x_mb_set.pop('miRNAseq')
    x_mb_set[2]=x_mb_set.pop('mRNAseq')
    x_mb_set[3]=x_mb_set.pop('RPPA')
    
    _, Lt, Lp, Lkl, Lps, Lkls, Lc = model.train(
        x_mb_set, y_mb, m_mb, alpha, beta, lr_rate, k_prob)

    tr_avg_Lt += Lt/STEPSIZE
    tr_avg_Lp += Lp/STEPSIZE
    tr_avg_Lkl += Lkl/STEPSIZE
    tr_avg_Lps += Lps/STEPSIZE
    tr_avg_Lkls += Lkls/STEPSIZE
    tr_avg_Lc += Lc/STEPSIZE
   
    x_mb_set, y_mb, m_mb = f_get_minibatch_set(
        min(np.shape(va_M)[0], mb_size), va_X_set, va_Y_onehot, va_M)
    
    x_mb_set[0]=x_mb_set.pop('Methylation')
    x_mb_set[1]=x_mb_set.pop('miRNAseq')
    x_mb_set[2]=x_mb_set.pop('mRNAseq')
    x_mb_set[3]=x_mb_set.pop('RPPA')
    
    Lt, Lp, Lkl, Lps, Lkls, Lc, _, _ = model.get_loss(
        x_mb_set, y_mb, m_mb, alpha, beta)

    va_avg_Lt += Lt/STEPSIZE
    va_avg_Lp += Lp/STEPSIZE
    va_avg_Lkl += Lkl/STEPSIZE
    va_avg_Lps += Lps/STEPSIZE
    va_avg_Lkls += Lkls/STEPSIZE
    va_avg_Lc += Lc/STEPSIZE 
    if (itr+1) % STEPSIZE == 0:
        
        va_X_set[0]=va_X_set.pop('Methylation')
        va_X_set[1]=va_X_set.pop('miRNAseq')
        va_X_set[2]=va_X_set.pop('mRNAseq')
        va_X_set[3]=va_X_set.pop('RPPA')
        y_pred, y_preds = model.predict_ys(va_X_set, va_M)
