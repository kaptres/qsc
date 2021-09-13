% author: mferhata
corr_matrix_c21 = [];
corr_matrix_c22 = [];

mydir = 'results/';
for file1=dir([mydir '*.txt'])'
    filename1       = [mydir file1.name];

    corr_with1_c21 = [];
    corr_with1_c22 = [];
    for file2=dir([mydir '*.txt'])'
        filename2           = [mydir file2.name];

        [c21s, c22s] = pairwise_correlation_collection2(filename1, filename2);
        corr_with1_c21 = [corr_with1_c21 c21s];
        corr_with1_c22 = [corr_with1_c22 c22s];
    end

    corr_matrix_c21 = [corr_matrix_c21; corr_with1_c21];
    corr_matrix_c22 = [corr_matrix_c22; corr_with1_c22];
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

for i=1:length(method_names); names{i} = char(method_names(i)); end;
names_rev = names(end:-1:1);
corr_matrix_c21 = corr_matrix_c21(end:-1:1,:);
corr_matrix_c22 = corr_matrix_c22(end:-1:1,:);

h = HeatMap (corr_matrix_c21, 'colormap', 'parula', 'rowlabels', names_rev, 'columnlabels', names, 'columnlabelsrotate', 45, 'DisplayRange', 2);
%h.addTitle ('collection2\_1');
h.plot; set (gcf, 'color', 'w', 'position', [1200 300 660 600]);
c = colorbar; set (c, 'Location', 'manual'); set (c, 'Position', [0.92 0.12 .03 .78]);
pause (1); 
children = get (gcf, 'children');
set (children(3), 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'XAxisLocation', 'top');
pause (.1); saveas (gcf, 'c21_corr.png');

h = HeatMap (corr_matrix_c22, 'colormap', 'parula', 'rowlabels', names_rev, 'columnlabels', names, 'columnlabelsrotate', 45, 'DisplayRange', 2);
%h.addTitle ('collection2\_2'); 
h.plot; set (gcf, 'color', 'w', 'position', [1200 300 660 600]);
c = colorbar; set (c, 'Location', 'manual'); set (c, 'Position', [0.92 0.12 .03 .78]);
pause (1);
children = get (gcf, 'children');
set (children(3), 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'XAxisLocation', 'top');
pause (.1); saveas (gcf, 'c22_corr.png'); 
