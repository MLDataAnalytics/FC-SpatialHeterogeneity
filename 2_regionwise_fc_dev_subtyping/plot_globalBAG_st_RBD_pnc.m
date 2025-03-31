clc;
clear;


% load analysis results
res_dir = '/cbica/home/lihon/comp_space_recover/ABCD_comp/dl_prediction/pred_res_cog_abcc_bl_cv/non_dl_pnc_res';
dat_dir = '/cbica/home/lihon/comp_space_recover/ABCD_comp/dl_prediction/data_cog_pnc_cv';

num_r = 400;

st_type = 'regional';
num_K = 3;

ro_covars = 1;
out_dir = [res_dir, filesep, 'vis_v2_r', num2str(num_r), '_ro', num2str(ro_covars)];

res_file = [out_dir, filesep, 'bag_res_K', num2str(num_K), '_', st_type, '_r', num2str(num_r), '.mat'];
load(res_file, 'st_idx', 'age_true', 'g_bag', 'c_cell');


% plot global BAG
newcolors = [0, 0.4470, 0.7410; ...
            0.4940 0.1840 0.5560; ...
            0.9290 0.6940 0.1250];

% RBD based subtype
h_st = figure;
colormap(newcolors);
hold on;
for ki=1:num_K
    k_idx = st_idx == ki;
    scatter(age_true(k_idx), g_bag(k_idx), 15, st_idx(k_idx), '+');
end
legend(c_cell, 'Location', 'best');
xlim([8, 23]);
xlabel('Chronological Age (Years)');
ylabel('Global BAG');

savefig(h_st, [out_dir, filesep, 'scatter_global_bag_all_K', num2str(num_K), '_', st_type, '_r', num2str(num_r), '.fig']);

disp('Finished.');

