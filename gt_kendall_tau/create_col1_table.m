head = 'Method & $\verb|w|=3$ & $\verb|w|=4$ & $\verb|w|=5$ & $\verb|c|=25$ & $\verb|c|=50$ & $\verb|c|=75$\\\midrule';
disp ('\begin{table}[h]');
disp ('\caption{\label{table:??}text}');
disp ('\begin{adjustbox}{width=\columnwidth,center}');
disp ('\centering');
disp ('\begin{tabular}{ccccccc}');
disp ('\toprule');
disp (head);
for i=1:length(cellarr)
    if (sum(cellarr{i}(:)) == 0) & (sum(cellarr2{i}(:)) == 0)
        continue;
    end
    string = sprintf('%35s&\t',names{i});
    for j=1:length(cellarr{1})
        if cellarr{i}(j)>=0
            score1 = sprintf('\\hphantom{-}%4.2f', cellarr{i}(j));
        else
            score1 = sprintf('%4.2f', cellarr{i}(j));
        end
        if cellarr2{i}(j)>=0
            score2 = sprintf('\\hphantom{-}%4.2f', cellarr2{i}(j));
        else
            score2 = sprintf('%4.2f', cellarr2{i}(j));
        end
        if j==length(cellarr{1})
            string = [string score1 '\hphantom{-}/' score2 '\\'];
        else
            string = [string score1 '\hphantom{-}/' score2 sprintf('&\t')];
        end
    end
    disp (string);
end
meanabs1 = mean ( abs(reshape([cellarr{:} ], [6 length(cellarr )])'));
meanabs2 = mean ( abs(reshape([cellarr2{:}], [6 length(cellarr2)])'));
meanabs_string = '\midrule $MA$ &';
for i=1:length(meanabs1);
    if meanabs1(i)>=0
        score1 = sprintf('\\hphantom{-}%4.2f', meanabs1(i));
    else
        score1 = sprintf('%4.2f', meanabs1(i));
    end
    if meanabs2(i)>=0
        score2 = sprintf('\\hphantom{-}%4.2f', meanabs2(i));
    else
        score2 = sprintf('%4.2f', meanabs2(i));
    end

    if i==length(meanabs1)
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
