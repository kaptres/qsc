% author: mferhata
% Es should be calculated before calling this (see calculate_E.m)
% t1 and t2 should be specified before calling
load_dirs = [...
        "~/DATA/SHREC_collection1/c+_Mar26" ...
        "~/DATA/SHREC_collection1/c-_Mar26" ...
        "~/DATA/SHREC_collection1/s+_Mar26" ...
        "~/DATA/SHREC_collection1/s-_Mar26" ...
        "~/DATA/SHREC_collection2/collection2_1_Mar26" ...
        "~/DATA/SHREC_collection2/collection2_2_Mar26" ...
        "~/DATA/SHREC_collection3/off_Mar26" ...
        ];
f = fopen (['results_' sprintf('%.2f',t1) '<t<=' sprintf('%.2f',t2) '.txt'], 'a');
bot_ind = floor (1 + t1 * 100);
top_ind = floor (1 + t2 * 100);
display(num2str(top_ind));

sss = {}; sss_index = 1;
for load_dir = load_dirs
    load_dir = char (load_dir);

    switch sss_index
    case 1
        fprintf (f, 'c+');
    case 2
        fprintf (f, 'c-');
    case 3
        fprintf (f, 's+');
    case 4
        fprintf (f, 's-');
    case 5
        fprintf (f, 'c21');
    case 6 
        fprintf (f, 'c22');
    case 7
        fprintf (f, 'off');
    end
    fprintf (f, '\n');

    ss = []; ss_index = 1;
    for file = dir (fullfile (load_dir, '*_Es.mat'))'
        load (fullfile (load_dir, file.name));

        ss(ss_index).name   = file.name;
        ss(ss_index).Es     = Es;
        s                   = mean(Es(bot_ind:top_ind));
        ss_index            = ss_index + 1;
        fprintf (f, '(%f) ', s);
    end
    sss{sss_index}  = ss;
    sss_index       = sss_index + 1;
    fprintf (f, '\n');
end
fclose (f);
