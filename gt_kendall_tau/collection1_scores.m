% author: mferhata
% computes the Kendall's tau corr. coeffs 
%   for the results of a method and the ground truth on Collection 1
function [cp_scores, cm_scores, sp_scores, sm_scores, cp_scores2, cm_scores2, sp_scores2, sm_scores2 ]...
        = collection1_scores (filename, disp_on)
    if ~exist ('disp_on', 'var')
        disp_on = false;
    end

    parse_submission2;

    if ~exist ('cp', 'var'); cp = zeros (450,1); end;
    if ~exist ('cm', 'var'); cm = zeros (450,1); end;
    if ~exist ('sp', 'var'); sp = zeros (450,1); end;
    if ~exist ('sm', 'var'); sm = zeros (450,1); end;
    if size(cp,2) == 2
        cp = cp(:,2);
        cm = cm(:,2);
        sp = sp(:,2);
        sm = sm(:,2);
    end

    cp  = reshape (cp, [9 50])';
    cm  = reshape (cm, [9 50])';
    sp  = reshape (sp, [9 50])';
    sm  = reshape (sm, [9 50])';

    cp_scores   = []; cp_scores2  = [];
    cm_scores   = []; cm_scores2  = [];
    sp_scores   = []; sp_scores2  = [];
    sm_scores   = []; sm_scores2  = [];

    nw  = [1:3; 4:6; 7:9];          % fix nws
    cnt = [1 4 7; 2 5 8; 3 6 9];    % fix counts
    qs  = [nw; cnt];
    for i=1:6
        cp_scores   = [cp_scores arrayfun(@(n) kendall_tau(cp(n, qs(i,:)), 1:3), 1:50)'];
        cm_scores   = [cm_scores arrayfun(@(n) kendall_tau(cm(n, qs(i,:)), 1:3), 1:50)'];
        sp_scores   = [sp_scores arrayfun(@(n) kendall_tau(sp(n, qs(i,:)), 1:3), 1:50)'];
        sm_scores   = [sm_scores arrayfun(@(n) kendall_tau(sm(n, qs(i,:)), 1:3), 1:50)'];
        % the following uses the MATLAB's built-in function for calculating the Kendall-tau
        %for n=1:50
        %    cp_scores2  = [cp_scores2; corr(cp(n, qs(i,:))', [1:3]', 'type', 'kendall')];
        %    cm_scores2  = [cm_scores2; corr(cm(n, qs(i,:))', [1:3]', 'type', 'kendall')];
        %    sp_scores2  = [sp_scores2; corr(sp(n, qs(i,:))', [1:3]', 'type', 'kendall')];
        %    sm_scores2  = [sm_scores2; corr(sm(n, qs(i,:))', [1:3]', 'type', 'kendall')];
        %end
    end
    %cp_scores2  = reshape (cp_scores2, [50 6]);
    %cm_scores2  = reshape (cm_scores2, [50 6]);
    %sp_scores2  = reshape (sp_scores2, [50 6]);
    %sm_scores2  = reshape (sm_scores2, [50 6]);

    [t, name, t] = fileparts(filename);
    if disp_on
        disp (name);
        disp ('       nw=3    nw=4    nw=5    cnt=25  cnt=50  cnt=75');
        disp (['c+:    ' sprintf('%.3f   ',sum(cp_scores/150))]);
        disp (['s+:    ' sprintf('%.3f   ',sum(sp_scores/150))]);
        %disp (['c-:    ' sprintf('%.3f   ',sum(cm_scores/150))]);
        %disp (['s-:    ' sprintf('%.3f   ',sum(sm_scores/150))]);
    end
end
