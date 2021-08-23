% author: mferhata
function [score, Es] = TIP_score (v)
    dt  = bwdist (~v, 'ch');
    dt  = dt / max(dt(:));
    Es  = [];
    edges = 0:1/1024:1;
    for d=0:100
        vals    = get_vals_from_ref (v, dt, d/100);
        N       = histcounts (vals, edges);
        N       = N / sum(N);
        E       = -sum (N(N>0) .* log2 (N(N>0)));
        Es      = [Es E];
    end
    score = mean (Es);
end
