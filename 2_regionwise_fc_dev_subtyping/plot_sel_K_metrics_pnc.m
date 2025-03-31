clc;
clear;


rt_dir = '/cbica/home/lihon/comp_space_recover/ABCD_comp/dl_prediction/pred_res_cog_abcc_bl_cv/non_dl_pnc_res/vis_v2_r400_ro1';
metric_mat_file = [rt_dir, filesep, 'sel_K_ver2_metric.mat'];

load(metric_mat_file);

m_ari = mean(ari_vec);
s_ari = std(ari_vec);

m_sumd = mean(md_vec);
s_sumd = std(md_vec);

h_f = figure; %hold on;
yyaxis left;
errorbar(k_set, m_ari, s_ari, '-o', 'MarkerFaceColor', '#0072BD');
xlabel('# of clusters');
ylabel('Reproducibility');

yyaxis right;
errorbar(k_set, m_sumd, s_sumd, '-s', 'MarkerFaceColor', '#D95319');
ylabel('Data fitting');

xlim([1.5, 7.5]);

fig_name = [rt_dir, filesep, 'sel_K_ver2_plot.fig'];
savefig(h_f, fig_name);

disp('Finished.');
