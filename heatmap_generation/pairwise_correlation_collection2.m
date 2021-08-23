% author: mferhata
% computes the pairwise correlations between two methods on collection1
function [c21s, c22s] = pairwise_correlation_collection2 (filename1, filename2)
    load_reverseds;
    %reverseds = [];

    filename = filename1;
    parse_submission2;
    [t, name, t] = fileparts (filename);
    if any (strcmp (name, reverseds));
        c21_1 = -c21; c22_1 = -c22;
    else
        c21_1 = c21; c22_1 = c22;
    end

    filename = filename2;
    parse_submission2;
    [t, name, t] = fileparts (filename);
    if any (strcmp (name, reverseds));
        c21_2 = -c21; c22_2 = -c22;
    else
        c21_2 = c21; c22_2 = c22;
    end

    c21s = corr(c21_1, c21_2, 'type', 'kendall');
    c22s = corr(c22_1, c22_2, 'type', 'kendall');
end
