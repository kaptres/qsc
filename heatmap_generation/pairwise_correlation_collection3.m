% author: mferhata
% computes the pairwise correlations between two methods on collection 3
% returns results acquired by   a) taking shapes in the same category into consideration
%                               b) taking all shapes into consideration
function [catdata, alldata] = pairwise_correlation_collection3 (filename1, filename2)
    load_reverseds;
    %reverseds = [];

    filename = filename1;
    parse_submission2;
    [t, name, t] = fileparts (filename);
    if any (strcmp (name, reverseds));
        off1 = -off;
    else
        off1 = off;
    end

    filename = filename2;
    parse_submission2;
    [t, name, t] = fileparts (filename);
    if any (strcmp (name, reverseds));
        off2 = -off;
    else
        off2 = off;
    end

    catdata = [];
    for i=1:19
        lo_ind = 1 + 20 * (i - 1);
        hi_ind = 20 * i;

        score = corr(off1(lo_ind:hi_ind), off2(lo_ind:hi_ind), 'type', 'kendall');
        %score = corr(off1(lo_ind:hi_ind), off2(lo_ind:hi_ind), 'type', 'spearman');
        catdata = [catdata score];
    end
    alldata = corr(off1, off2, 'type', 'kendall');
    %alldata = corr(off1, off2, 'type', 'spearman');
end
