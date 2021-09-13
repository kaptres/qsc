% author: mferhata
% computes the Kendall's tau corr. coeffs 
%   for the results of a method and the ground truth on Collection 2
function [c21s, c22s, results1, results2] = collection2_scores (filename, disp_on)
    if ~exist ('disp_on', 'var')
        disp_on = false;
    end

    parse_submission2;

    if ~exist ('c21', 'var')
        if disp_on; display (filename); end;
        return
    end

    [t,name,t] = fileparts (filename);
    if disp_on
        disp (name);
        disp(['  |--c21 ' sprintf('G%d  \t',1:5) 'Sum']);
    end

    clear G;
    for k=1:1%size(c21,2)
        temp = c21(:,k);
        G{1}= temp([16 14 12 17 18 19 20]);
        G{2}= temp([15 13 24 22]);
        G{3}= temp([23 21 25]);
        G{4}= temp([11 10 7 1 2]);
        G{5}= temp([9 8 4 5 6 3]);

        results = [];
        for i=1:5; 
            s       = kendall_tau (G{i}, 1:length(G{i}));
            s       = s / (length(G{i}) * (length(G{i})-1) / 2);
            results = [results s];
        end
        c21s    = sum(results);
        results = [results c21s];
        results1= results;


        if disp_on
            disp (['  |      ' sprintf('%.2f\t', results)]);
        end
    end

    if disp_on
        disp (['  |--c22 ' sprintf('G%d  \t',1:6) 'Sum']);
    end
    clear G;
    for k=1:1%size(c22,2)
        temp = c22(:,k);
        temp = temp - min(temp(:));
        temp = temp / max(temp(:));
        G{1}= temp([18 24 17 20]);
        G{2}= temp([25 21 22 19]);
        G{3}= temp([5 3 2]);
        G{4}= temp([23 16 15]);
        G{5}= temp([6 4 14 12 13]);
        G{6}= temp([7 9 10 1 11 8]);

        results = [];
        for i=1:6; 
            s       = 0;
            for j=1:length(G{i})
                for k=j+1:length(G{i});
                    s = s + abs (G{i}(j) - G{i}(k));
                end
            end
            results = [results s];
        end
        c22s    = sum(results);
        results = [results c22s];
        results2= results;


        if disp_on
            disp (['  |      ' sprintf('%.2f\t', results)]);
        end
    end
end
