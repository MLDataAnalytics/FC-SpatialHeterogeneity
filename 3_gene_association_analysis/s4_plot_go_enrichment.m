clc;
clear;


enrich_file = 'C:\Users\lihon\Desktop\abcd_study\pnc\vis_v2_r400_ro1\GOCOMPONENT.xls.xlsx';
tbl = readtable(enrich_file, 'VariableNamingRule', 'preserve');

go_term = tbl{:, 1};
go_desc = tbl{:, 2};
p_val = tbl{:, 3};
q_val = tbl{:, 4};
enrichment = tbl{:, 5};
B = tbl{:, 7};
b = tbl{:, 9};

figure; hold on;
x_val = enrichment;
y_val = -log10(q_val);
sz = b / 3; 
s = scatter(x_val, y_val, sz, 'filled', 'MarkerEdgecolor', 'k');
s.AlphaData = -log10(q_val);
s.MarkerFaceAlpha = 'flat';
xlabel('Enrichment');
ylabel('-log_1_0(q-val)');

for i=1:min(10, length(x_val))
        yrand = rand(1) * 0.02 * ((-1).^(i-1));
        text(x_val(i)*1.01, y_val(i)*(1+yrand), go_desc{i}, 'FontSize', 7);
end

disp('Finished.');


