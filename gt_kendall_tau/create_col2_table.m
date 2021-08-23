head = 'Method & Group 1 & Group 2 & Group 3 & Group 4 & Group 5 & Group 6 & Sum\\\midrule';
disp ('\begin{table}[h]');
disp ('\caption{\label{table:??}text}');
disp ('\begin{adjustbox}{width=\columnwidth,center}');
disp ('\centering');
disp ('\begin{tabular}{cccccccc}');
disp ('\toprule');
disp (head);
for i=1:length(G1s)
    string = sprintf('%35s&\t',names{i});
    for j=1:length(G2s{1})
        if j==length(G1s{1})
            score1 = sprintf('\\hphantom{-}--');
        elseif j==length(G2s{1})
            score1 = normalize_num (G1s{i}(j-1),j, G2s);
        else
            score1 = normalize_num (G1s{i}(j),j, G2s);
        end
        score2 = normalize_num (G2s{i}(j),j, G2s);
        if j==length(G2s{1})
            string = [string score1 '\hphantom{-}/' score2 '\\'];
        else
            string = [string score1 '\hphantom{-}/' score2 sprintf('&\t')];
        end
    end
    if i==length(cellarr)
        string = [string '\bottomrule'];
    end
    disp (string);
end
meanabs1 = mean ( abs(reshape([G1s{:}], [6 length(G1s)])'));
meanabs2 = mean ( abs(reshape([G2s{:}], [7 length(G2s)])'));
meanabs_string = '\midrule $MA$ &';
for i=1:length(meanabs2);
    if i==length(meanabs1)
        score1 = '--';
    elseif i==length(meanabs1)+1
        score1 = sprintf('\\hphantom{-}%4.2f', meanabs1(i-1));
    elseif meanabs1(i)>=0
        score1 = sprintf('\\hphantom{-}%4.2f', meanabs1(i));
    else
        score1 = sprintf('%4.2f', meanabs1(i));
    end
    if meanabs2(i)>=0
        score2 = sprintf('\\hphantom{-}%4.2f', meanabs2(i));
    else
        score2 = sprintf('%4.2f', meanabs2(i));
    end

    if i==length(meanabs2)
        meanabs_string = [meanabs_string score1 '\hphantom{-}/' score2 '\\'];
    else
        meanabs_string = [meanabs_string score1 '\hphantom{-}/' score2 sprintf('&\t')];
    end
end
meanabs_string = [meanabs_string '\bottomrule'];

disp (meanabs_string);
disp ('\end{tabular}');
disp ('\end{adjustbox}');
disp ('\end{table}');

function string = normalize_num (n, j, G2s)
    if j == length(G2s{1})
        if floor(abs(n)/10) == 0
            if n>=0
                string = sprintf('\\hphantom{1}\\hphantom{-}%4.2f',n);
            else
                string = sprintf('\\hphantom{-}%4.2f',n);
            end
        else
            if n>=0
                string = sprintf('\\hphantom{-}%4.2f',n);
            else
                string = sprintf('%4.2f',n);
            end

        end
    else
        if n>=0
            string = sprintf('\\hphantom{-}%4.2f',n);
        else
            string = sprintf('%4.2f',n);
        end
    end
end
