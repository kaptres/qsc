% author: mferhata
% computes Kendall's tau
function res = kendall_tau (p, q)
    res = 0;
    for i=1:length(p)-1
        for j=i+1:length(p)
            if ( (p(i) - p(j)) * (q(i) - q(j)) ) > 0
                res = res + 1;
            elseif ( (p(i) - p(j)) * (q(i) - q(j)) ) < 0
                res = res - 1;
            end
        end
    end
end
