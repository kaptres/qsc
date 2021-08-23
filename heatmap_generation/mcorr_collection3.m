% author: mferhata
corr_matrix_catdata = [];
corr_matrix_alldata = [];

mydir = 'results/';
for file1=dir([mydir '*.txt'])'
    filename1       = [mydir file1.name];

    corr_with1_catdata = [];
    corr_with1_alldata = [];
    for file2=dir([mydir '*.txt'])'
        filename2           = [mydir file2.name];

        [catdata, alldata]  = pairwise_correlation_collection3(filename1, filename2);
        corr_with1_catdata  = [corr_with1_catdata mean(catdata)];
        corr_with1_alldata  = [corr_with1_alldata alldata];
    end

    corr_matrix_catdata = [corr_matrix_catdata; corr_with1_catdata];
    corr_matrix_alldata = [corr_matrix_alldata; corr_with1_alldata];
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
corr_matrix_catdata = corr_matrix_catdata(end:-1:1,:);
corr_matrix_alldata = corr_matrix_alldata(end:-1:1,:);

h = HeatMap (corr_matrix_catdata, 'colormap', 'parula', 'rowlabels', names_rev, 'columnlabels', names, 'columnlabelsrotate', 45, 'displayrange', 2);
%h.addTitle ('collection3: means of category-wise correlations'); 
h.plot; set (gcf, 'color', 'w', 'position', [600 300 660 600]);
c = colorbar; set (c, 'Location', 'manual'); set (c, 'Position', [0.92 0.12 .03 .78]);
pause (1); 
children = get (gcf, 'children'); 
set (children(3), 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'XAxisLocation', 'top');
pause (.1); export_fig c3_cat_corr.pdf; 

h = HeatMap (corr_matrix_alldata, 'colormap', 'parula', 'rowlabels', names_rev, 'columnlabels', names, 'columnlabelsrotate', 45, 'DisplayRange', 2);
%h.addTitle ('collection3: correlations over all shapes');
h.plot; set (gcf, 'color', 'w', 'position', [600 300 660 600]);
c = colorbar; set (c, 'Location', 'manual'); set (c, 'Position', [0.92 0.12 .03 .78]);
pause (1); 
children = get (gcf, 'children'); 
set (children(3), 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'XAxisLocation', 'top');
pause (.1); export_fig c3_all_corr.pdf; 
