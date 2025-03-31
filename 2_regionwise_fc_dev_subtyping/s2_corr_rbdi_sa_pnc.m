clc;
clear;


atlas = cifti_read('C:\Users\lihon\Desktop\abcd_study\atlas\Schaefer\Schaefer2018_400Parcels_17Networks_order.dlabel.nii');
atlas_parcel = atlas.cdata;
atlas_label = unique(atlas_parcel(atlas_parcel>0));
num_parcel = length(atlas_label);

% load regional brain development (RBD) index patterns
K = 3;

cen_dir = 'C:\Users\lihon\Desktop\abcd_study\pnc\vis_v2_r400_ro1';
rbag_cen = cell(K, 1);
for ki=1:K
        rbag_cen{ki} = cifti_read([cen_dir, filesep, 'regional_bag_st_cen_', num2str(ki), 'in', num2str(K), '.dscalar.nii']);
end

vert_ind_l = rbag_cen{1}.diminfo{1}.models{1}.vertlist + 1;
vert_ind_r = rbag_cen{1}.diminfo{1}.models{2}.vertlist + 1 + 32492;
vert_ind = [vert_ind_l'; vert_ind_r'];
atlas_parcel_v2 = atlas.cdata(vert_ind,:);

rbag_mat = zeros(num_parcel, K);
for ai=1:num_parcel
        ai_ind = atlas_parcel_v2==atlas_label(ai);
        
        for ki=1:K
                rbag_mat(ai, ki) = mean(rbag_cen{ki}.cdata(ai_ind));
        end
end

% load S-A axis rank data 
sa_file = 'C:\Users\lihon\Desktop\abcd_study\S-A_ArchetypalAxis-main\Schaefer400_17Network\Sensorimotor_Association_Axis_AverageRanks.csv';
sa_mat = readmatrix(sa_file);

sa_vec = sa_mat(:,2);   % final rank
corr_sa_rbag = corr(sa_vec, rbag_mat);

% spin test
spin_test = load('C:\Users\lihon\Desktop\abcd_study\pnc\spin_perm_parcellation\spin_perm_id_S400_N17.mat');

num_rot = size(spin_test.perm_id, 2);
rot_cor = zeros(num_rot, size(rbag_mat,2));
for ri=1:num_rot
        rot_ind = spin_test.perm_id(:, ri);
        rot_sa = sa_vec(rot_ind);
        
        rot_cor(ri, :) = corr(rot_sa, rbag_mat);
end

spin_p = zeros(size(corr_sa_rbag));
for bi=1:length(corr_sa_rbag)
        bi_real = corr_sa_rbag(bi);
        bi_rot = rot_cor(:, bi);
        if bi_real > 0
                bi_p = sum(bi_rot>bi_real) / num_rot;
        else
                bi_p = sum(bi_rot<bi_real) / num_rot;
        end
        spin_p(bi) = bi_p;
end

%% gen plot
for ki=1:K
        figure; hold on; 
        rbag_vec = rbag_mat(:, ki);
        scatter(sa_vec, rbag_vec, 20, [248, 222, 126]./255, 'filled');
        Fit = polyfit(sa_vec, rbag_vec, 1);
        plot(sa_vec, polyval(Fit, sa_vec), '-', 'Color', '#00A36C', 'LineWidth', 1.5); 

        xlabel('S-A Axis Rank');
        ylabel('Regional BAG');
        xlim([0, 400]);
        
        y_lim_val = [-1, 0; -1, 1; -0.1, 1.2];
        ylim(y_lim_val(ki,:));
        if spin_p(ki)==0
                spin_str = ' \itp<0.0001';
        elseif spin_p(ki)>=0.001
                spin_str = [' \itp=', sprintf('%0.3f', spin_p(ki))];
        else
                spin_str = [' \itp=', num2str(spin_p(ki))];
        end
        
        y_pos = [-0.1, -0.85, 0];
        text(50, y_pos(ki), ['\itr=', sprintf('%.3f', corr_sa_rbag(ki)), spin_str]);
        
        %figure; hold on;
        h_pos = [0.65, 0.2, 0.25, 0.25; ...
                         0.65, 0.2, 0.25, 0.25; ...
                         0.65, 0.65, 0.25, 0.25];
        axes('pos', h_pos(ki,:));
        histogram(rot_cor(:,ki), 20, 'FaceColor', [229, 228, 226]./255);
        set(gca, 'ytick', []);
        xline(corr_sa_rbag(ki), 'Color', '#00A36C', 'LineWidth', 1.5);
        box off;
end
disp('Finished.');

