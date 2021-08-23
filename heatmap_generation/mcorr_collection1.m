% author: mferhata
corr_matrix_cp = [];
corr_matrix_cm = [];
corr_matrix_sp = [];
corr_matrix_sm = [];

mydir = 'results/';
for file1=dir([mydir '*.txt'])'
    corr_with1_cp   = [];
    corr_with1_cm   = [];
    corr_with1_sp   = [];
    corr_with1_sm   = [];

    filename1       = [mydir file1.name];
    for file2=dir([mydir '*.txt'])'
        filename2           = [mydir file2.name];
        [cps, cms, sps, sms]= pairwise_correlation_collection1 (filename1, filename2);

        corr_with1_cp   = [corr_with1_cp mean(cps)];
        corr_with1_cm   = [corr_with1_cm mean(cms)];
        corr_with1_sp   = [corr_with1_sp mean(sps)];
        corr_with1_sm   = [corr_with1_sm mean(sms)];
    end

    corr_matrix_cp = [corr_matrix_cp; corr_with1_cp];
    corr_matrix_cm = [corr_matrix_cm; corr_with1_cm];
    corr_matrix_sp = [corr_matrix_sp; corr_with1_sp];
    corr_matrix_sm = [corr_matrix_sm; corr_with1_sm];
end

method_names = [ ...
            "[16]" "[14]" "[21]-1" "[21]-2" ...
            "[15]" ...
            "[29]" "[28]" ...
            "[27]" "$\mathcal{C}_{CRE}$" "[4]" ...
            "[26]" "$\mathcal{C}_{\sigma}$" "$\mathcal{C}_1$" "$\mathcal{C}_{PC}$" ...
            "[24]" "[23]" "$\mathcal{C}_2$" "[22]-1" "[22]-2" ...
            ];
%method_names = method_names([1 2 3 4 13 5 9 12 8 17 18 19 16 6 11 15 7 10 14]);

for i=1:length(method_names); names{i} = char(method_names(i)); end
names_rev = names(end:-1:1);
corr_matrix_cp = corr_matrix_cp(end:-1:1,:);
corr_matrix_sp = corr_matrix_sp(end:-1:1,:);
corr_matrix_cm = corr_matrix_cm(4:-1:1,1:4);
corr_matrix_sm = corr_matrix_sm(4:-1:1,1:4);

h = HeatMap (corr_matrix_cp, 'colormap', 'parula', 'rowlabels', names_rev, 'columnlabels', names, 'columnlabelsrotate', 45, 'DisplayRange', 2);
%h.addTitle ('c+');
h.plot; set (gcf, 'color', 'w', 'position', [600 300 660 600]);
c = colorbar; set (c, 'Location', 'manual'); set (c, 'Position', [0.92 0.12 .03 .78]);
pause (1); 
children = get (gcf, 'children'); 
set (children(3), 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'XAxisLocation', 'top');
pause (1); saveas (gcf, 'cp_corr.png'); 

h = HeatMap (corr_matrix_sp, 'colormap', 'parula', 'rowlabels', names_rev, 'columnlabels', names, 'columnlabelsrotate', 45, 'DisplayRange', 2);
%h.addTitle ('s+');
h.plot; set (gcf, 'color', 'w', 'position', [600 300 660 600]);
c = colorbar; set (c, 'Location', 'manual'); set (c, 'Position', [0.92 0.12 .03 .78]);
pause (1); 
children = get (gcf, 'children'); 
set (children(3), 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'XAxisLocation', 'top');
pause (1); saveas (gcf, 'sp_corr.png'); 

h = HeatMap (corr_matrix_cm, 'colormap', 'parula', 'rowlabels', names(4:-1:1), 'columnlabels', names(1:4), 'columnlabelsrotate', 0, 'DisplayRange', 2);
%h.addTitle ('c-');
h.plot; set (gcf, 'color', 'w', 'position', [600 300 660 600]);
c = colorbar; set (c, 'Location', 'manual'); set (c, 'Position', [0.92 0.12 .03 .78]);
pause (1); 
children = get (gcf, 'children'); 
set (children(3), 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'XAxisLocation', 'top');
pause (1); saveas (gcf, 'cm_corr.png'); 

h = HeatMap (corr_matrix_sm, 'colormap', 'parula', 'rowlabels', names(4:-1:1), 'columnlabels', names(1:4), 'columnlabelsrotate',0, 'DisplayRange', 2);
%h.addTitle ('s-');
h.plot; set (gcf, 'color', 'w', 'position', [600 300 660 600]);
c = colorbar; set (c, 'Location', 'manual'); set (c, 'Position', [0.92 0.12 .03 .78]);
pause (1); 
children = get (gcf, 'children'); 
set (children(3), 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'XAxisLocation', 'top');
pause (1); saveas (gcf, 'sm_corr.png'); 
