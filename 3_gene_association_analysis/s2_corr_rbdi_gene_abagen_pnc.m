clc;
clear;

addpath('C:\work\code\download\freesurfer\matlab');
addpath('C:\work\code\download\cifti-matlab-master');


% load brain atlas
num_roi = 400;
num_roi_str = num2str(num_roi);

atlas_file = ['C:\Users\lihon\Desktop\abcd_study\atlas\Schaefer2018_', num_roi_str, 'Parcels_17Networks_order.dlabel.nii'];
atlas_cii = cifti_read(atlas_file);
atlas_label = unique(atlas_cii.cdata(atlas_cii.cdata>0));
num_parcel = length(atlas_label);

% load gene expression data from abagen
gene_file = ['C:\Users\lihon\Desktop\AHBA\abagen\abagen_gene_Schaefer2018_', num_roi_str, 'Parcels_17Networks_stable_gene.csv'];
gene_tbl = readtable(gene_file, 'VariableNamingRule', 'preserve', 'ReadRowNames', true);

gene_symbol = gene_tbl.Properties.VariableNames;
gene_dat = single(table2array(gene_tbl));
gene_dat = gene_dat(1:num_roi/2, :);    % left hemi only
clear gene_tbl;

% keep brain-express genes (from Seidlitz et. al, Transcriptomic and cellular decoding of regional brain vulnerability to neurogenetic disorders)
brain_gene_file = 'C:\Users\lihon\Desktop\abcd_study\pnc\brain_genes_HPA.txt';
bgfid = fopen(brain_gene_file, 'r');
brain_gene_lst = textscan(bgfid, '%s');
brain_gene_lst = brain_gene_lst{1};
fclose(bgfid);

sel_be_gene = zeros(length(gene_symbol), 1);
for si=1:length(gene_symbol)
        sind = find(strcmp(gene_symbol{si}, brain_gene_lst));
        if ~isempty(sind)
                sel_be_gene(si) = 1;
        end
end
sel_be_gene = sel_be_gene > 0;

gene_dat = gene_dat(:, sel_be_gene);
gene_symbol = gene_symbol(sel_be_gene);

% remove nan gene and brain region
sel_roi = sum(isnan(gene_dat), 2) < size(gene_dat,2);
gene_dat = gene_dat(sel_roi, :);
sel_g = sum(isnan(gene_dat)) == 0;
gene_dat = gene_dat(:, sel_g);

gene_dat = zscore(gene_dat);
gene_symbol = gene_symbol(sel_g);


% load regional brain development (RBD) index patterns
K = 3;
cen_dir = ['C:\Users\lihon\Desktop\abcd_study\pnc\vis_v2_r', num_roi_str, '_ro1'];
rbag_cen = cell(K, 1);
for ki=1:K
        rbag_cen{ki} = cifti_read([cen_dir, filesep, 'regional_bag_st_cen_', num2str(ki), 'in', num2str(K), '.dscalar.nii']);
end

vert_ind_l = rbag_cen{1}.diminfo{1}.models{1}.vertlist + 1;
vert_ind_r = rbag_cen{1}.diminfo{1}.models{2}.vertlist + 1 + 32492;
vert_ind = [vert_ind_l'; vert_ind_r'];
atlas_parcel = atlas_cii.cdata(vert_ind,:);

rbag_mat = zeros(num_parcel, K);
for ai=1:num_parcel
        ai_ind = atlas_parcel==atlas_label(ai);
        
        for ki=1:K
                rbag_mat(ai, ki) = mean(rbag_cen{ki}.cdata(ai_ind));
        end
end

% pattern P2, left hemisphere
lh_c2 = rbag_mat(1:num_roi/2, 2);
lh_c2 = lh_c2(sel_roi, :);

% corr ana
[gene_corr, gene_p] = corr(gene_dat, lh_c2);

% ranked list
[sort_gcor_val, sort_gcor_ind] = sort(gene_corr, 'descend');
sort_gene_symbol = gene_symbol(sort_gcor_ind);

% sig genes by spin-test
spin_test = load(['C:\Users\lihon\Desktop\abcd_study\pnc\spin_perm_parcellation\spin_perm_id_S', num_roi_str, '_N17.mat']);
perm_id = spin_test.perm_id(1:num_roi/2, :);

num_rot = size(perm_id, 2);
rot_cor = zeros(num_rot, size(gene_dat,2));
lh_c2_real = rbag_mat(1:num_roi/2,2);
for ri=1:num_rot
        rot_ind = perm_id(:, ri);
        rot_lh_c2 = lh_c2_real(rot_ind);
        rot_lh_c2 = rot_lh_c2(sel_roi,:);
        
        rot_cor(ri, :) = corr(gene_dat, rot_lh_c2);
end

spin_p = zeros(size(gene_corr));
for bi=1:length(gene_corr)
        bi_real = gene_corr(bi);
        bi_rot = rot_cor(:, bi);
        if bi_real > 0
                bi_p = sum(bi_rot>bi_real) / num_rot;
        else
                bi_p = sum(bi_rot<bi_real) / num_rot;
        end
        spin_p(bi) = bi_p;
end
spin_p_fdr = mafdr(spin_p, 'BHFDR', true);

%% save csv for cell type enrichment
% AHBAdata file to get entrez_ID, gene_index
ahba_mat = load("C:\Users\lihon\Desktop\abcd_study\pnc\gene_decoding\sex_differences_adolescence\Data\gene_decoding\AHBAdata.mat");

entrezId1 = zeros(size(sort_gene_symbol)) - 1;
for gi=1:length(sort_gene_symbol)
        gind = find(strcmp(sort_gene_symbol{gi}, ahba_mat.probeInformation.GeneSymbol));
        if ~isempty(gind)
                entrezId1(gi) = ahba_mat.probeInformation.EntrezID(gind);
        else
                disp([num2str(gi), '. ', sort_gene_symbol{gi}]);
        end
end
g_avail = entrezId1 > 0;
gene_sb1 = sort_gene_symbol(g_avail);
entrezId1 = entrezId1(g_avail);
z1 = sort_gcor_val(g_avail);

tab_all = table(gene_sb1', categorical(entrezId1)', z1);
tab_all.Properties.VariableNames = {'gene_name', 'entrez_ID', 'z_uncorr'};

out_file = [cen_dir, filesep, 'corr_rbag_gene_pnc.csv'];
writetable(tab_all, out_file);

%% save sig gene list
p_thr = 0.05;
p_thr_str = '5e-2';
sig_gene = spin_p_fdr < p_thr;
sig_gene_symbol = gene_symbol(sig_gene);
sig_gene_cor = gene_corr(sig_gene);

figure;
ct_mat = [0.8500 0.3250 0.0980; 0 0.4470 0.7410];
c_vec = ones(size(sig_gene_cor));
c_vec(sig_gene_cor<=0) = 2;
c_mat = ct_mat(c_vec, :);
wordcloud(sig_gene_symbol, abs(sig_gene_cor), 'color', c_mat);

% save gene list
out_file = [cen_dir, filesep, 'sig_gene_lst_cen2_p', p_thr_str, '_spin.txt'];
out_file_pos = [cen_dir, filesep, 'sig_gene_lst_cen2_p', p_thr_str, '_spin_pos.txt'];
out_file_neg =  [cen_dir, filesep, 'sig_gene_lst_cen2_p', p_thr_str, '_spin_neg.txt'];
ofid = fopen(out_file, 'w');
ofid_p = fopen(out_file_pos, 'w');
ofid_n = fopen(out_file_neg, 'w');
for gi=1:length(sig_gene_symbol)
        fprintf(ofid, '%s\n', sig_gene_symbol{gi});
        if sig_gene_cor(gi)>0
                fprintf(ofid_p, '%s\n', sig_gene_symbol{gi});
        else
                fprintf(ofid_n, '%s\n', sig_gene_symbol{gi});
        end
end
fclose(ofid);
fclose(ofid_p);
fclose(ofid_n);

ref_file = [cen_dir, filesep, 'ref_gene_lst.txt'];
ofid = fopen(ref_file, 'w');
for gi=1:length(gene_symbol)
        fprintf(ofid, '%s\n', gene_symbol{gi});
end
fclose(ofid);

disp('Finished.');

