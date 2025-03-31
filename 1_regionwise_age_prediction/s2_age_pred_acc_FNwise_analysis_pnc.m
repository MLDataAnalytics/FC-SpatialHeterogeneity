clc;
clear;


addpath('/cbica/home/lihon/code/Download/cifti-matlab-master');


res_dir = '/cbica/home/lihon/comp_space_recover/ABCD_comp/dl_prediction/pred_res_cog_abcc_bl_cv/non_dl_pnc_res';
dat_dir = '/cbica/home/lihon/comp_space_recover/ABCD_comp/dl_prediction/data_cog_pnc_cv';

num_r = 400;

ro_covars = 1;
out_dir = [res_dir, filesep, 'vis_v2_r', num2str(num_r), '_ro', num2str(ro_covars)];

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

num_fd = 5;

age_true = [];
age_pred = [];
sex = [];
motion = [];
for fi=1:num_fd
    fi_res = load([res_dir, filesep, 'pred_res_age_fold', num2str(fi-1), '_Ridge_rw.mat']);
    
    age_true = [age_true; fi_res.y_true'];
    age_pred = [age_pred; fi_res.y_pred_mat];
    
    fi_dat = load([dat_dir, filesep, 'pnc_fc_fold', num2str(fi-1), '.mat']);
    sex = [sex; fi_dat.tes_dat(:,end-1)];
    motion = [motion; fi_dat.tes_dat(:,end)];
end

%% regional performance
cor_vec = partialcorr(age_true, age_pred, [sex, motion]); 
mae_vec = mean(abs(repmat(age_true,1, size(age_pred,2))-age_pred));
r_acc_mat = [cor_vec; mae_vec];
r_acc_cell = {'partialcor', 'mae'};

% save metric map to cifti
tpl_lab_file = ['/cbica/home/lihon/atlas/Schaefer2018/Schaefer2018_', num2str(num_r), 'Parcels_17Networks_order.dlabel.nii'];
tpl_lab = cifti_read(tpl_lab_file);

tpl_scalar_file = '/cbica/home/lihon/comp_space_recover/hcp_mn/hcp_d_mn/more_measures/HCD0001305_V1_MR/HCD0001305_V1_MR.ArealDistortion_sm10_MSMAll.32k_fs_LR.dscalar.nii';
tpl_scalar = cifti_read(tpl_scalar_file);

for ai = 1:size(r_acc_mat,1)
    ai_vec = r_acc_mat(ai,:);
    att_cdata = zeros(size(tpl_lab.cdata));
    for pi=1:length(ai_vec)
        att_cdata(tpl_lab.cdata==pi) = ai_vec(pi);
    end
    m_att_cifti = tpl_scalar;
    att_l = att_cdata(1:32492,:);
    att_l = att_l(tpl_scalar.diminfo{1}.models{1}.vertlist+1,:);
    att_r = att_cdata(32492+1:64984,:);
    att_r = att_r(tpl_scalar.diminfo{1}.models{2}.vertlist+1,:);
    m_att_cifti.cdata = [att_l; att_r];
    m_att_file = [out_dir, filesep, 'regional_', r_acc_cell{ai}, '.dscalar.nii'];
    cifti_write(m_att_cifti, m_att_file);
end

%% get region to FN info
[num, txt, raw] = xlsread('/cbica/home/lihon/atlas/Schaefer2018/atlas-Schaefer2018v0143_desc-400ParcelsAllNetworks_dseg.xlsx');

index_17 = num(:,5);
network_label_7 = txt(2:end,3);

[label_set, ia, ic] = unique(network_label_7);
label_set_copy = label_set;
label_set{1} = 'Frontoparietal';
label_set{3} = 'Dorsal Attention';
label_set{5} = 'Ventral Attention';
label_set{6} = 'Somatomotor';
label_set{7} = 'Visual';

lab_color = [230, 148, 34; ...
             205, 62, 78; ...
             0, 118, 14; ...
             220, 248, 164; ...
             196, 58, 250; ...
             70, 130, 180; ...
             120, 18, 134] ./ 255;

% boxchart
fn_idx = ic;
fn_idx = categorical(fn_idx, 1:length(label_set), label_set);

cor_vec_ = cor_vec(index_17);
mae_vec_ = mae_vec(index_17);

h_cor = figure; 
hold on;
for fi=1:length(label_set)
    boxchart(fn_idx(ic==fi), cor_vec_(ic==fi), 'BoxWidth', 0.3, 'BoxFaceColor', lab_color(fi,:)); 
end
colororder(lab_color); 
ylabel('Accuracy (r)');
savefig(h_cor, [out_dir, filesep, 'regional_partial_cor_fnwise.fig']);

h_mae = figure; 
hold on;
for fi=1:length(label_set)
    boxchart(fn_idx(ic==fi), mae_vec_(ic==fi), 'BoxWidth', 0.3, 'BoxFaceColor', lab_color(fi,:)); 
end
colororder(lab_color);
ylabel('MAE (years)');
savefig(h_mae, [out_dir, filesep, 'regional_mae_fnwise.fig']);

disp('Finished.');


