% author: mferhata
% computes the Kendall's tau corr. coeffs 
%   for the results of a method and the ground truth on Collection 3
function [cat_mean, all_mean, cat_std, all_std] = collection3_scores (filename, disp_on)
    if ~exist ('disp_on', 'var')
        disp_on = false;
    end
    parse_submission2;
    segmentations;
    mean_scores     = reshape ( arrayfun (@(n) mean(human_segmentations {n}   ), 1:length(human_segmentations)), [20 19]);
    std_scores      = reshape ( arrayfun (@(n) std (human_segmentations {n}   ), 1:length(human_segmentations)), [20 19]);
    mad_scores0     = reshape ( arrayfun (@(n) mad (human_segmentations {n}, 0), 1:length(human_segmentations)), [20 19]);
    mad_scores1     = reshape ( arrayfun (@(n) mad (human_segmentations {n}, 1), 1:length(human_segmentations)), [20 19]);
    %mean_scores_2   = reshape ( arrayfun (@(n) mean(human_segmentations2{n}   ), 1:length(human_segmentations)), [20 19]);
    %std_scores_2    = reshape ( arrayfun (@(n) std (human_segmentations2{n}   ), 1:length(human_segmentations)), [20 19]);
    %mad_scores0_2   = reshape ( arrayfun (@(n) mad (human_segmentations2{n}, 0), 1:length(human_segmentations)), [20 19]);
    %mad_scores1_2   = reshape ( arrayfun (@(n) mad (human_segmentations2{n}, 1), 1:length(human_segmentations)), [20 19]);
    all_scores      = {mean_scores  std_scores};%  mad_scores0  mad_scores1  mean_scores_2  std_scores_2  mad_scores0_2  mad_scores1_2};

    [t, name, t] = fileparts (filename);
    if disp_on; disp (name); end;
    for j=1:1
        if j==1
            if disp_on; disp ('    |All segmentations included'); end;
            scores  = all_scores(1:2);%4);
        else
            if disp_on; disp ('    |Excluding some segmentations'); end;
            scores  = all_scores(5:8);
        end
        for k=1:2
            if k==1
                if disp_on; disp ('    |----Mean'); end;
            elseif k==2
                if disp_on; disp ('    |----STD'); end;
            elseif k==3
                if disp_on; disp ('    |----Mean absolute deviation'); end;
            elseif k==4
                if disp_on; disp ('    |----Median absolute deviation'); end;
            end
            score = scores{k};
            a = [];
            for i=1:19
                a = [a corr(off(20*i-19:20*i,1), score(:,i), 'type', 'kendall')];
            end
            if disp_on; disp('    |   |----category no:          1         2       3      4       5        6      7        8      9       10       11    12        13     14       15      16     17       18     19'); end;

            if disp_on; disp (['    |   |----for each category:   ' sprintf('%.3f\t', a)]); end;
            allscore = corr (off(:,1), score(:), 'type', 'kendall');
            %allscore = allscore * 2 / (379 * 380);
            if disp_on; disp (['    |   |----over all dataset:    '  sprintf('%.3f', allscore)]); end;


            if k==1; cat_mean = mean(a); all_mean = allscore; end
            if k==2; cat_std  = mean(a); all_std  = allscore; end
        end
    end
end
