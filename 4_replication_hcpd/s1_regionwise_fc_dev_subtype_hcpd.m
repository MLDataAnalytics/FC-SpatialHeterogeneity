clc;
clear;

addpath('/cbica/home/lihon/code/Download/cifti-matlab-master');


num_r = 400;

mat_dir = '/cbica/home/lihon/comp_space_recover/ABCD_comp/dl_prediction/data_cog_hcpd_cv_combined';
res_dir = ['/cbica/home/lihon/comp_space_recover/ABCD_comp/dl_prediction/pred_res_cog_abcc_bl_cv/non_dl_hcpd_res_r', num2str(num_r), '_combined_xcpd'];

ro_covars = 1;
out_dir = [res_dir, filesep, 'vis_v2_r', num2str(num_r), '_ro', num2str(ro_covars), '_updated']; %with site

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

num_fd = 5;

age_true = [];
age_pred = [];
sex = [];
site = [];
motion = [];
for fi=1:num_fd
    fi_res = load([res_dir, filesep, 'pred_res_age_fold', num2str(fi-1), '_Ridge_rw.mat']);
    
    age_true = [age_true; fi_res.y_true'];
    age_pred = [age_pred; fi_res.y_pred_mat];
    
    cog_dat = load([mat_dir, filesep, 'hcpd_fc_fold', num2str(fi-1), '.mat']);
    sex = [sex; cog_dat.tes_dat(:,end-2)];
    site = [site; cog_dat.tes_dat(:,end-1)];
    motion = [motion; cog_dat.tes_dat(:,end)];
end
site_c = categorical(site);
site_oh = onehotencode(site_c, 2);
site_oh = site_oh(:,1:end-1);

%% regional performance
cor_vec = corr(age_true, age_pred);
mae_vec = mean(abs(repmat(age_true,1, size(age_pred,2))-age_pred));
r_acc_mat = [cor_vec; mae_vec];
r_acc_cell = {'cor', 'mae'};

tpl_lab_file = '/cbica/home/lihon/atlas/Schaefer2018/Schaefer2018_400Parcels_17Networks_order.dlabel.nii';
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

%% regional development index
num_roi = size(age_pred, 2);
bag_mat = zeros(size(age_pred));
for ri=1:num_roi
    ri_y = age_pred(:, ri);
    
    % bias correction, regress out covariates
    ri_b = glmfit([age_true, motion, sex, site_oh], ri_y);
    ri_bag = ri_y - ri_b(2).*age_true - ri_b(3).*motion - ri_b(4).*sex - site_oh*ri_b(5:7) - ri_b(1);
    
    bag_mat(:,ri) = ri_bag;
end
    
%% subtyping
% subtyping using PNC RBD index patterns
st_type = 'regional';
num_K = 3;

pnc_cen_dir = '/cbica/home/lihon/comp_space_recover/ABCD_comp/dl_prediction/pred_res_cog_abcc_bl_cv/non_dl_pnc_res/vis_v2_r400_ro1';
pnc_cen_file = [pnc_cen_dir, '/bag_cen_K', num2str(num_K), '_', st_type, '_r', num2str(num_r), '.mat'];
pnc_cen = load(pnc_cen_file);

st_cen = pnc_cen.st_cen;

[~, st_idx] = pdist2(st_cen, bag_mat, 'cosine', 'Smallest', 1);
st_idx = st_idx';

c_cell = cell(num_K, 1);
for ki=1:num_K
    c_cell{ki} = ['P', num2str(ki)];
end

% eval subtypes in terms of cognitive measures
f1 = [];
f2 = [];
f3 = [];
for fi=1:num_fd
    cog_dat = load([mat_dir, filesep, 'hcpd_fc_fold', num2str(fi-1), '.mat']);

    f1 = [f1; cog_dat.tes_dat(:,5)];    % fluid cog
    f2 = [f2; cog_dat.tes_dat(:,6)];    % crystal cog
    f3 = [f3; cog_dat.tes_dat(:,8)];    % total cog
end
cog_fid = [ones(size(f1)); 2*ones(size(f2)); 3*ones(size(f3))];
cog_fid = categorical(cog_fid, 1:3, {'Fluid cognition', 'Crystal cognition', 'Total cognition'});
cog_f = [f1; f2; f3];
grp_idx = [st_idx; st_idx; st_idx];
grp_idx = categorical(grp_idx, 1:num_K, c_cell);

h_bc = figure; boxchart(cog_fid, cog_f, 'GroupByColor', grp_idx); ylabel('Score');
legend(c_cell, 'Location', 'northoutside', 'NumColumns', 3); legend('boxoff');
hold on; for ci=2:3; h_x = xline(ci-0.5, 'Color', [0.5,0.5,0.5], 'HandleVisibility','off'); end
savefig(h_bc, [out_dir, filesep, 'boxchart_bag_cog_', st_type, '_K', num2str(num_K), '.fig']);

cog_mat = [f1, f2, f3];
num_cog = size(cog_mat, 2);
h_mat = ones(num_K,num_K,num_cog);
p_mat = ones(num_K,num_K,num_cog);
for fi=1:num_cog
    for ci=1:num_K
        for cj=ci+1:num_K
            [p_mat(ci,cj,fi), h_mat(ci,cj,fi)] = ranksum(cog_mat(st_idx==ci,fi), cog_mat(st_idx==cj,fi));
        end
    end
end

% test sex difference across subtypes
[tbl, chi2, s_p] = crosstab(st_idx, sex);

% test age, motion difference across sbutypes
m_h = ones(num_K);
m_p = ones(num_K);
a_h = ones(num_K);
a_p = ones(num_K);
for ci=1:num_K
    for cj=ci+1:num_K
        [m_p(ci,cj), m_h(ci,cj)] = ranksum(motion(st_idx==ci,:), motion(st_idx==cj,:));
        [a_p(ci,cj), a_h(ci,cj)] = ranksum(age_true(st_idx==ci,:), age_true(st_idx==cj,:));
    end
end

res_file = [out_dir, filesep, 'bag_res_K', num2str(num_K), '_', st_type, '_r', num2str(num_r), '.mat'];
save(res_file);

disp('Finished.');

