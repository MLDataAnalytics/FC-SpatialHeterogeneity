clc;
clear;

addpath('C:\work\code\download\cifti-matlab-master');


bag_res_file = 'C:\Users\lihon\Desktop\abcd_study\pnc\vis_v2_r400_ro1\bag_cen_K3_regional_r400.mat';
bag_res = load(bag_res_file);

bag_p2 = bag_res.bag_mat(bag_res.st_idx==2, :);
bag_p3 = bag_res.bag_mat(bag_res.st_idx==3, :);

num_r = size(bag_p2, 2);

test_23 = zeros(num_r, 2);
test_32 = zeros(num_r, 2);
test_two = zeros(num_r, 2);

for ri=1:num_r
        [test_two(ri, 1), h_, stat_] = ranksum(bag_p2(:, ri), bag_p3(:, ri));
        
        [test_23(ri, 1), h_23, stat_23] = ranksum(bag_p2(:, ri), bag_p3(:, ri), 'tail', 'right');
        [test_32(ri, 1), h_32, stat_32] = ranksum(bag_p2(:, ri), bag_p3(:, ri), 'tail', 'left');
        
        test_23(ri, 2) = stat_23.zval;
        test_32(ri, 2) = stat_32.zval;
end
val_23 = mafdr(test_23(:, 1), 'BHFDR', true);
val_32 = mafdr(test_32(:, 1), 'BHFDR', true);

val_two = mafdr(test_two(:, 1), 'BHFDR', true);
val_two = single(val_two<0.05);

% save res to cifti
tpl_lab_file = 'C:\Users\lihon\Desktop\abcd_study\atlas\Schaefer2018_400Parcels_17Networks_order.dlabel.nii';
tpl_lab = cifti_read(tpl_lab_file);

tpl_scalar_file = 'C:\Users\lihon\Desktop\abcd_study\pnc\vis_v2_r400_ro1\regional_cor.dscalar.nii';
tpl_scalar = cifti_read(tpl_scalar_file);

val_mat = [];
val_mat(1, :) = test_23(:, 2);
val_mat(2, :) = single(val_23<0.05);
val_mat(3, :) = single(val_32<0.05);
val_mat(4, :) = val_mat(2, :) + val_mat(3, :)*-1;

fdr_str = {'z', 'p2gt3', 'p3gt2', 'com23'};
out_dir = 'C:\Users\lihon\Desktop\abcd_study\pnc\vis_v2_r400_ro1';
save([out_dir, filesep, 'bag_compare_p2_vs_p3.mat'], 'val_mat', 'val_two');

for ai = 1:size(val_mat,1)
    ai_vec = val_mat(ai,:);
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
    m_att_file = [out_dir, filesep, 'bag_compare_val_', fdr_str{ai}, '.dscalar.nii'];
    cifti_write(m_att_cifti, m_att_file);
end

disp('Finished.');
