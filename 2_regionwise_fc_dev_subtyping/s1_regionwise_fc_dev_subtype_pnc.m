clc;
clear;

addpath('/cbica/home/lihon/code/Download/cifti-matlab-master');


res_dir = '/cbica/home/lihon/comp_space_recover/ABCD_comp/dl_prediction/pred_res_cog_abcc_bl_cv/non_dl_pnc_res';
dat_dir = '/cbica/home/lihon/comp_space_recover/ABCD_comp/dl_prediction/data_cog_pnc_cv';

num_r = 400;

ro_covars = 1;  % 0: bag; 1: bag r/o sex, motion; 2: bag r/o sex, motion, g-bag
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

%% regionwise prediction performance
cor_vec = corr(age_true, age_pred);
mae_vec = mean(abs(repmat(age_true,1, size(age_pred,2))-age_pred));
r_acc_mat = [cor_vec; mae_vec];
r_acc_cell = {'cor', 'mae'};

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

%% global and regional brain development index
% global development index, predicted brain age from whole FC
m_y = [];
for fi=1:num_fd
    fi_res = load([res_dir, filesep, 'pred_res_age_fold', num2str(fi-1), '_Ridge.mat']);
    m_y = [m_y; fi_res.y_pred'];
end

% bias correction, regress out covariates
m_b = glmfit([age_true, motion, sex], m_y);
g_bag = m_y - m_b(2).*age_true - m_b(3).*motion - m_b(4).*sex - m_b(1);

% regional development index
num_roi = size(age_pred, 2);
bag_mat = zeros(size(age_pred));
bag_mat_orig = bag_mat;
for ri=1:num_roi
    ri_y = age_pred(:, ri);
    
    % bias correction, regress out covariates
    ri_b = glmfit([age_true, motion, sex], ri_y);
    ri_bag = ri_y - ri_b(2).*age_true - ri_b(3).*motion - ri_b(4).*sex - ri_b(1);

    bag_mat_orig(:,ri) = ri_bag;
    bag_mat(:,ri) = ri_bag;
end

%% select K (number of subtypes)
sel_K = 0;

if sel_K > 0
    rng(17);
    
    k_set = 2:7;
    num_rep = 100;

    ari_vec = zeros(num_rep, length(k_set));
    md_vec = zeros(num_rep, length(k_set));
    for rpi=1:num_rep
        disp(num2str(rpi));

        num_sbj = size(bag_mat, 1);
        rpi_ind = randperm(num_sbj);
        half_ind = round(num_sbj/2);

        bag_mat1 = bag_mat(rpi_ind(1:half_ind),:);
        bag_mat2 = bag_mat(rpi_ind(half_ind+1:end),:);

        for ki=1:length(k_set)
            [ki_idx1, ki_cen1, ki_sumd1] = kmeans(bag_mat1, k_set(ki), 'Distance', 'cosine', 'Replicates', 100);
            [ki_idx2, ki_cen2, ki_sumd2] = kmeans(bag_mat2, k_set(ki), 'Distance', 'cosine', 'Replicates', 100);

            [~, ki_idx1_2] = pdist2(ki_cen2, bag_mat1, 'cosine', 'Smallest', 1);
            [~, ki_idx2_1] = pdist2(ki_cen1, bag_mat2, 'cosine', 'Smallest', 1);

            for dki=1:k_set(ki)
                ki_sumd1(dki) = ki_sumd1(dki) / sum(ki_idx1==dki);
                ki_sumd2(dki) = ki_sumd2(dki) / sum(ki_idx2==dki);
            end
            md_vec(rpi, ki) = 1 - (mean(ki_sumd1) + mean(ki_sumd2)) / 2;

            % ari
            ki_ari_1 = RandIndex(ki_idx1, ki_idx1_2);
            ki_ari_2 = RandIndex(ki_idx2, ki_idx2_1);
            ari_vec(rpi, ki) = (ki_ari_1 + ki_ari_2) / 2;
        end
    end
    save([out_dir, filesep, 'sel_K_metric.mat'], 'k_set', 'ari_vec', 'md_vec');

    k_cell = cell(1, length(k_set));
    for ki=1:length(k_set)
        k_cell{ki} = num2str(k_set(ki));
    end
    h_sk1 = figure; boxplot(ari_vec, 'Labels', k_cell); xlabel('# of cluster'); ylabel('ARI');
    savefig(h_sk1, [out_dir, filesep, 'boxplot_sel_K_ARI.fig']);
    h_sk2 = figure; boxplot(md_vec, 'Labels', k_cell); xlabel('# of cluster'); ylabel('Within-cluster dist');
    savefig(h_sk2, [out_dir, filesep, 'boxplot_sel_K_sumd.fig']);
end

%% subtyping
st_type = 'regional'; %'regional' or 'global'
num_K = 3;
num_rep = 100;

rng('default');
rng(17);

if strcmp(st_type, 'regional')
    [st_idx, st_cen] = kmeans(bag_mat, num_K, 'Distance', 'cosine', 'Replicates', num_rep);
    st_cen_orig = zeros(size(st_cen));
    for ki=1:num_K
        st_cen_orig(ki,:) = mean(bag_mat_orig(st_idx==ki,:));
    end

    c_cell = cell(num_K, 1);
    for ki=1:num_K
        c_cell{ki} = ['P', num2str(ki)];
    end
elseif strcmp(st_type, 'global')
    [st_idx, st_cen] = kmeans(g_bag, num_K, 'Replicates', num_rep);
    st_cen_orig = st_cen;
    
    c_cell = cell(num_K, 1);
    for ki=1:num_K
        c_cell{ki} = ['P', num2str(ki), '-G'];
    end
else
    error(['unrecognized st_type: ', st_type]);
end

% sort idx by g_bag
st_g_bag = zeros(num_K,1);
for ki=1:num_K
    st_g_bag(ki) = mean(g_bag(st_idx==ki));
end
[~, sort_ind] = sort(st_g_bag);
st_cen = st_cen(sort_ind,:);
st_cen_orig = st_cen_orig(sort_ind,:);
st_idx_org = st_idx;
for si=1:length(sort_ind)
    st_idx(st_idx_org==sort_ind(si)) = si;
end

cen_file = [out_dir, filesep, 'bag_cen_K', num2str(num_K), '_', st_type, '_r', num2str(num_r), '.mat'];
save(cen_file, 'st_cen', 'st_idx', 'bag_mat', 'st_cen_orig', 'bag_mat_orig', 'g_bag');

if strcmp(st_type, 'regional')
    for ki=1:num_K
        ki_cen = st_cen_orig(ki,:);
        
        att_cdata = zeros(size(tpl_lab.cdata));
        for pi=1:length(ki_cen)
            att_cdata(tpl_lab.cdata==pi) = ki_cen(pi);
        end

        m_att_cifti = tpl_scalar;
        att_l = att_cdata(1:32492,:);
        att_l = att_l(tpl_scalar.diminfo{1}.models{1}.vertlist+1,:);
        att_r = att_cdata(32492+1:64984,:);
        att_r = att_r(tpl_scalar.diminfo{1}.models{2}.vertlist+1,:);
        m_att_cifti.cdata = [att_l; att_r];
        m_att_file = [out_dir, filesep, st_type, '_bag_st_cen_', num2str(ki), 'in', num2str(num_K), '.dscalar.nii'];
        cifti_write(m_att_cifti, m_att_file);
    end
end

% eval subtypes in terms of cognitive measures
f1 = [];
f2 = [];
f3 = [];
for fi=1:num_fd
    fi_res = load([res_dir, filesep, 'pred_res_c1_fold', num2str(fi-1), '_Ridge_rw.mat']);
    f1 = [f1; fi_res.y_true'];
    
    fi_res = load([res_dir, filesep, 'pred_res_c2_fold', num2str(fi-1), '_Ridge_rw.mat']);
    f2 = [f2; fi_res.y_true'];
    
    fi_res = load([res_dir, filesep, 'pred_res_c3_fold', num2str(fi-1), '_Ridge_rw.mat']);
    f3 = [f3; fi_res.y_true'];
end
cog_fid = [ones(size(f1)); 2*ones(size(f2)); 3*ones(size(f3))];
cog_fid = categorical(cog_fid, 1:3, {'Executive', 'Social', 'Memory'});    % Executive cognition, social cognition, episodic memory 
cog_f = [f1; f2; f3];
grp_idx = [st_idx; st_idx; st_idx];
grp_idx = categorical(grp_idx, 1:num_K, c_cell);

h_bc = figure; boxchart(cog_fid, cog_f, 'GroupByColor', grp_idx); ylabel('Score'); ylim([-4, 4]);
legend(c_cell, 'Location', 'northoutside', 'NumColumns', 3); legend('boxoff');
hold on; for ci=2:3; h_x = xline(ci-0.5, 'Color', [0.5,0.5,0.5], 'HandleVisibility','off'); end
savefig(h_bc, [out_dir, filesep, 'boxchart_bag_cog_', st_type, '_K', num2str(num_K), '.fig']);

cog_mat = [f1, f2, f3];
h_mat = ones(num_K,num_K,3);
p_mat = ones(num_K,num_K,3);

for fi=1:3
    for ci=1:num_K
        for cj=ci+1:num_K
            ci_dat = cog_mat(st_idx==ci,fi);
            cj_dat = cog_mat(st_idx==cj,fi);
            
            if ~isempty(ci_dat) && ~isempty(cj_dat)
                [p_mat(ci,cj,fi), h_mat(ci,cj,fi)] = ranksum(ci_dat, cj_dat);
            end
        end
    end
end

% test sex difference across subtypes
[tbl, chi2, s_p] = crosstab(st_idx, sex);

% test age, g_bag, motion difference across sbutypes
m_h = ones(num_K);
m_p = ones(num_K);
a_h = ones(num_K);
a_p = ones(num_K);
bag_h = ones(num_K);
bag_p = ones(num_K);
for ci=1:num_K
    for cj=ci+1:num_K
        [m_p(ci,cj), m_h(ci,cj)] = ranksum(motion(st_idx==ci,:), motion(st_idx==cj,:));
        [a_p(ci,cj), a_h(ci,cj)] = ranksum(age_true(st_idx==ci,:), age_true(st_idx==cj,:));
        [bag_p(ci,cj), bag_h(ci,cj)] = ranksum(g_bag(st_idx==ci,:), g_bag(st_idx==cj,:));
    end
end

res_file = [out_dir, filesep, 'bag_res_K', num2str(num_K), '_', st_type, '_r', num2str(num_r), '.mat'];
save(res_file);

disp('Finished.');

