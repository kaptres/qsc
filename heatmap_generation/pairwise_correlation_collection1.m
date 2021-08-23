% author: mferhata
% computes the pairwise correlations between two methods on collection1
function [cps, cms, sps, sms] = pairwise_correlation_collection1 (filename1, filename2)
    load_reverseds;
    %reverseds = [];

    filename = filename1;
    parse_submission2;
    if ~exist ('cp', 'var'); cp1 = zeros (450,1); else cp1 = cp; end;
    if ~exist ('cm', 'var'); cm1 = zeros (450,1); else cm1 = cm; end;
    if ~exist ('sp', 'var'); sp1 = zeros (450,1); else sp1 = sp; end;
    if ~exist ('sm', 'var'); sm1 = zeros (450,1); else sm1 = sm; end;
    [t, name, t] = fileparts (filename);
    if any (strcmp (name, reverseds))
        cp1 = -cp1; cm1 = -cm1; sp1 = -sp1; sm1 = -sm1;
    end

    filename = filename2;
    parse_submission2;
    if ~exist ('cp', 'var'); cp2 = zeros (450,1); else cp2 = cp; end;
    if ~exist ('cm', 'var'); cm2 = zeros (450,1); else cm2 = cm; end;
    if ~exist ('sp', 'var'); sp2 = zeros (450,1); else sp2 = sp; end;
    if ~exist ('sm', 'var'); sm2 = zeros (450,1); else sm2 = sm; end;
    [t, name, t] = fileparts (filename);
    if any (strcmp (name, reverseds))
        cp2 = -cp2; cm2 = -cm2; sp2 = -sp2; sm2 = -sm2;
    end

    cps = []; cms = []; sps = []; sms = [];
    for i=1:50
        lo_i = 1 + 9 * (i - 1);
        hi_i = 9 * i;
        cps = [cps; corr(cp1(lo_i:hi_i), cp2(lo_i:hi_i), 'type', 'kendall')];
        cms = [cms; corr(cm1(lo_i:hi_i), cm2(lo_i:hi_i), 'type', 'kendall')];
        sps = [sps; corr(sp1(lo_i:hi_i), sp2(lo_i:hi_i), 'type', 'kendall')];
        sms = [sms; corr(sm1(lo_i:hi_i), sm2(lo_i:hi_i), 'type', 'kendall')];
    end
end
