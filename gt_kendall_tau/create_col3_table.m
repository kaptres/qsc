head = 'Method & $\mu_{\text{cat}}$ & $\mu_{\text{all}}$ & $\sigma_{\text{cat}}$ & $\sigma_{\text{all}}$\\\midrule';
disp ('\begin{table}[h]');
disp ('\caption{\label{table:??}text}');
disp ('\centering');
disp ('\begin{tabular}{ccccc}');
disp ('\toprule');
disp (head);
scores = reshape([col3s{:}], [4 length(col3s)])';
for i=1:length(col3s)
    string = sprintf('%35s&\t',names{i});
    for j=1:4
        if col3s{i}(j) >= 0
            score = sprintf ('\\hphantom{-}%4.3f\t', col3s{i}(j));
        else
            score = sprintf ('%4.3f\t', col3s{i}(j));
        end
        if j==4
            string = [string score '\\'];
        else
            string = [string score '&'];
        end
    end
    disp(string);
end
meanabs_string = sprintf('%35s&\t','\midrule$MA$');
meanabs = mean (abs (scores));
for j=1:4
    if meanabs(j) >= 0
        score = sprintf ('\\hphantom{-}%4.3f\t', meanabs(j));
    else
        score = sprintf ('%4.3f\t', meanabs(j));
    end
    if j==4
        meanabs_string = [meanabs_string score '\\\bottomrule'];
    else
        meanabs_string = [meanabs_string score '&'];
    end
end
disp (meanabs_string);
disp ('\end{tabular}');
disp ('\end{table}');

function string = normalize_num (n, j, G2s)
    if j == length(G2s{1})
        if floor(abs(n)/10) == 0
            if n>=0
                string = sprintf('\\hphantom{1}\\hphantom{-}%4.3f',n);
            else
                string = sprintf('\\hphantom{-}%4.3f',n);
            end
        else
            if n>=0
                string = sprintf('\\hphantom{-}%4.3f',n);
            else
                string = sprintf('%4.3f',n);
            end

        end
    else
        if n>=0
            string = sprintf('\\hphantom{-}%4.3f',n);
        else
            string = sprintf('%4.3f',n);
        end
    end
end
