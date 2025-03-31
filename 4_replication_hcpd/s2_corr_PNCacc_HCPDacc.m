clc;
clear;


atlas = cifti_read('C:\Users\lihon\Desktop\abcd_study\atlas\Schaefer\Schaefer2018_400Parcels_17Networks_order.dlabel.nii');
atlas_parcel = atlas.cdata;
atlas_label = unique(atlas_parcel(atlas_parcel>0));
num_parcel = length(atlas_label);

% load regionwise prediction accuracy on PNC and HCP-D
pnc_dir = 'C:\Users\lihon\Desktop\abcd_study\pnc\vis_v2_r400_ro1';
pnc_acc = cifti_read([pnc_dir, filesep, 'regional_cor.dscalar.nii']);

hcpd_dir = 'C:\Users\lihon\Desktop\abcd_study\pnc\hcpd_replicate\non_dl_hcpd_res_r400_combined_xcpd\vis_v2_r400_ro1_updated';
hcpd_acc = cifti_read([hcpd_dir, filesep, 'regional_cor.dscalar.nii']);

vert_ind_l = pnc_acc.diminfo{1}.models{1}.vertlist + 1;
vert_ind_r = pnc_acc.diminfo{1}.models{2}.vertlist + 1 + 32492;
vert_ind = [vert_ind_l'; vert_ind_r'];
atlas_parcel_v2 = atlas.cdata(vert_ind,:);

pnc_acc_vec = zeros(num_parcel, 1);
hcpd_acc_vec = zeros(num_parcel, 1);
for ai=1:num_parcel
        ai_ind = atlas_parcel_v2==atlas_label(ai);
        
        pnc_acc_vec(ai) = mean(pnc_acc.cdata(ai_ind));
        hcpd_acc_vec(ai) = mean(hcpd_acc.cdata(ai_ind));
end

corr_pnc_hcpd = corr(pnc_acc_vec, hcpd_acc_vec);

% spin test
spin_test = load('C:\Users\lihon\Desktop\abcd_study\pnc\spin_perm_parcellation\spin_perm_id_S400_N17.mat');

num_rot = size(spin_test.perm_id, 2);
rot_cor = zeros(num_rot, 1);
for ri=1:num_rot
        rot_ind = spin_test.perm_id(:, ri);
        rot_acc = hcpd_acc_vec(rot_ind);
        
        rot_cor(ri, :) = corr(rot_acc, pnc_acc_vec);
end

bi_real = corr_pnc_hcpd;
bi_rot = rot_cor;
if bi_real > 0
        bi_p = sum(bi_rot>bi_real) / num_rot;
else
        bi_p = sum(bi_rot<bi_real) / num_rot;
end
spin_p = bi_p;

%% gen plot
figure; hold on; 
scatter(pnc_acc_vec, hcpd_acc_vec, 20, [253, 218, 13]./255, 'filled');
Fit = polyfit(pnc_acc_vec, hcpd_acc_vec, 1);
acc_range = min(pnc_acc_vec):0.001:max(pnc_acc_vec);
plot(acc_range, polyval(Fit, acc_range), '-', 'Color',  '#00A36C', 'LineWidth', 1.5); 
xlabel('PNC Acc (\itr)');
ylabel('HCP-D Acc (\itr)');
xlim([0.2, 0.6]);
ylim([0.2, 0.6]);
if spin_p==0
        spin_str = ' \itp<0.0001';
else
        spin_str = [' \itp=', num2str(spin_p)];
end
text(0.225, 0.225, ['\itr=', sprintf('%.3f', corr_pnc_hcpd), spin_str]);

%figure; hold on;
axes('pos', [0.175, 0.65, 0.25, 0.25]);
histogram(rot_cor, 20, 'FaceColor', [229, 228, 226]./255);
set(gca, 'ytick', []);
xline(corr_pnc_hcpd, 'Color',  '#00A36C', 'LineWidth', 1.5);
box off;

disp('Finished.');

