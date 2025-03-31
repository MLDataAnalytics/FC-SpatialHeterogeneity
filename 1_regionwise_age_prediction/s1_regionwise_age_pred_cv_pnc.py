import os
import numpy as np

import pandas as pd

from sklearn.linear_model import RidgeCV
import scipy.io as sio


splits_dir = r"/cbica/home/lihon/comp_space/ABCD_comp/dl_prediction/data_cog_pnc_cv"
# 5-fold cross-validation
data_id = ["fold0", "fold1", "fold2", "fold3", "fold4"]
data_dir = r"/cbica/home/lihon/comp_space/bbl_pnc_resting/pnc_parcel_data_r400"

res_dir = r"/cbica/home/lihon/comp_space/ABCD_comp/dl_prediction/pred_res_cog_abcc_bl_cv/non_dl_pnc_res"
if not os.path.isdir(res_dir):
    os.makedirs(res_dir)

# region-wise FC profile
fea_type = "fc_p2p"     

# age prediction
y_ind = 4
y_type = "age"          

# regularization parameter in ridge regression
alphas_set = np.exp2(np.arange(16) - 10)

# number of cortical regions
num_r = 400

for fi in np.arange(len(data_id)):
    print("Fold " + str(fi) + " ...")
    tra_file = splits_dir + "/pnc_fc_tra_" + data_id[fi] + ".csv"
    tes_file = splits_dir + "/pnc_fc_tes_" + data_id[fi] + ".csv"

    # load data
    tra_dat_labels = pd.read_csv(tra_file, header=None)
    tes_dat_labels = pd.read_csv(tes_file, header=None)

    print("   load training data ...")
    num_tra = len(tra_dat_labels)
    x_tra = np.zeros((num_tra, num_r, num_r))
    y_tra = np.zeros((num_tra,))
    for di in np.arange(num_tra):
        di_path = os.path.join(data_dir, tra_dat_labels.iloc[di,0])

        di_data = sio.loadmat(di_path)[fea_type]
        di_data[np.isnan(di_data)] = 0.0

        di_label = tra_dat_labels.iloc[di, y_ind]

        x_tra[di] = di_data
        y_tra[di] = di_label

    print("   load testing data ...")
    num_tes = len(tes_dat_labels)
    x_tes = np.zeros((num_tes, num_r, num_r))
    y_tes = np.zeros((num_tes,))
    for di in np.arange(num_tes):
        di_path = os.path.join(data_dir, tes_dat_labels.iloc[di,0])

        di_data = sio.loadmat(di_path)[fea_type]
        di_data[np.isnan(di_data)] = 0.0

        di_label = tes_dat_labels.iloc[di, y_ind]

        x_tes[di] = di_data
        y_tes[di] = di_label

    # ridge regression   
    print("   model training ...")
    pred_y_mat = np.zeros((num_tes, num_r))

    for ri in np.arange(num_r):
        clf_p2p = RidgeCV(alphas=alphas_set).fit(x_tra[:,ri,:], y_tra)
        ri_pred_y_tes = clf_p2p.predict(x_tes[:,ri,:])
        pred_y_mat[:,ri] = ri_pred_y_tes
    pred_y_tes = np.mean(pred_y_mat, 1)
    
    fi_mae = np.mean(np.abs(y_tes-pred_y_tes))
    fi_corr = np.corrcoef(y_tes, pred_y_tes)[0,1]
    print("  p2p: MAE " + str(fi_mae) + " corr " + str(fi_corr))
    
    fi_out_file = res_dir + r"/pred_res_" + y_type + "_" + data_id[fi] + "_Ridge_rw.mat"
    
    sio.savemat(fi_out_file, {"y_true":y_tes, "y_pred":pred_y_tes, "y_pred_mat":pred_y_mat,
                              "mae":fi_mae, "corr":fi_corr})

print("Finished.")

