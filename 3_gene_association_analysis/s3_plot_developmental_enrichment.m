clc;
clear;


txt = 'C:\Users\lihon\Desktop\abcd_study\pnc\vis_v2_r400_ro1\SEA_brain_region_and_development_sigGene5e-2.txt';
tbl = readtable(txt);

struc_stage = tbl{:, 1};
p_val = tbl{:, 2};
p_val_fdr = tbl{:, 3};

struc_all = cell(size(struc_stage));
stage_all = cell(size(struc_stage));
for si=1:length(struc_stage)
        si_cell = strsplit(struc_stage{si}, '.');
        struc_all{si} = si_cell{1};
        stage_all{si} = strjoin(si_cell(2:end), ' ');
end

[struc_set, ia1, ic1] = unique(struc_all);
[stage_set, ia2, ic2] = unique(stage_all);

num_struc = length(struc_set);
num_stage = length(stage_set);
p_mat = zeros(num_struc, num_stage);
for i=1:length(ic1)
        p_mat(ic1(i), ic2(i)) = p_val_fdr(i);
end
log_p_mat = log10(p_mat);

stage_order = [9, 7, 1, 2, 4, 6, 3, 8, 5, 10];
[~, sort_stage_ind] = sort(stage_order);
stage_sort = stage_set(sort_stage_ind);
log_p_mat_sort = -log_p_mat(:, sort_stage_ind);

figure; 
imagesc(log_p_mat_sort, [0, 6]); %colormap('cool');
c = colorbar; c.Label.String = '-log_1_0(q-val)';
xticks(1:num_stage);
xticklabels(stage_sort);
yticks(1:num_struc);
yticklabels(struc_set);

% plot for cortex only
ctx_val = log_p_mat_sort(3, :);
ctx_val_late = ctx_val(7:end);
stage_late = stage_sort(7:end);

figure;
bar(ctx_val_late);
xticks(1:length(ctx_val_late));
xticklabels(stage_late);
ylabel('-log_1_0(q-val)');
ylim([0, 7]);
box off;

disp('Finished.');

