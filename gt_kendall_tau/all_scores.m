% author: mferhata
% creates all of the results by calling the scripts collection[1,2,3]_scores ()
% the resulting cell-variable names contains hardcoded latex references
names = {};
col11 = {};
col12 = {};
col13 = {};
col14 = {};
col21 = {};
col22 = {};
col31 = {};
col32 = {};
G1s   = {};
G2s   = {};
col3s = {};

i = 1;
for file=dir('results/*.txt')'
    filename = ['results/' file.name];

    disp(filename);
    [cp, cm, sp, sm]        = collection1_scores (filename);
    [c21, c22, G1, G2]      = collection2_scores (filename);
    [c31, c32, c33, c34]    = collection3_scores (filename);
    
    % Collection 1 scores of the methods
    col11{i} = sum(cp)/150; % scores for additively perturbed cubes
    col12{i} = sum(cm)/150; % scores for subtractively perturbed cubes
    col13{i} = sum(sp)/150; % scores for additively perturbed spheres
    col14{i} = sum(sm)/150; % scores for subtractively perturbed spheres

    % Collection 2 scores of the methods
    G1s{i}   = G1;          % scores for the first family
    G2s{i}   = G2;          % scores for the second family

    % Collection 3 scores of the methods
    % (two GTs: mean and std, category-wise and across all shapes)
    col3s{i} = [c31 c32 c33 c34];


    [t, name, t] = fileparts (file.name);
    names {i} = map_filenames_to_citations (name);
    i = i + 1;
end

function citation = map_filenames_to_citations (name)
    citation = '\cite{';
    switch name
    case "01_Ferhat"
        citation = [citation 'arslan2020'];
    case "02_MuratGenctav"
        citation = [citation 'Genctav2016'];
    case "03_Gardiner_Brassey_1"
        citation = [citation 'Gardiner2018'];
    case "04_Gardiner_Brassey_2"
        citation = [citation 'Gardiner2018'];
    case "05_AsliGenctav"
        citation = [citation 'genctavtari2019'];
    case "08_fractal_dim_massRobust"
        citation = [citation 'borkowski1999fractal'];
    case "09_fractal_dim2"
        citation = [citation 'hayward1989three'];
    case "10_k_regularity"
        citation = [citation 'VasselleG93'];
    case "13_complexity1"
        citation = '{${\mathcal C}_{CRE}$';
    case "13_complexityPage"
        citation = [citation 'page2003shape'];
    case "14_complexity5NormG16"
        citation = [citation 'matsumoto2019quantification'];
    case "15_complexity3NormG16"
        citation = '{${\mathcal C}_{\sigma}$';
    case "16_convexA"
        citation = '{$\mathcal C_1$';
    case "17_complexityView1"
        citation = '{${\mathcal C}_{PC}$';
    case "18_complexityBrinkhoff"
        citation = [citation 'brinkhoff1995measuring'];
    case "18_convexJ"
        citation = [citation 'zunic'];
    case "19_convexL"
        citation = '{$\mathcal C_2$';
    case "20_convexificationArea"
        citation = [citation 'rosin-lesions'];
    case "21_convexificationRevArea"
        citation = [citation 'rosin-lesions'];
    end
    citation = [citation '}'];
    switch name
    case "03_Gardiner_Brassey_1"
        citation = [citation '-1'];
    case "04_Gardiner_Brassey_2"
        citation = [citation '-2'];
    case "20_convexificationArea"
        citation = [citation '-1'];
    case "21_convexificationRevArea"
        citation = [citation '-2'];
    end
    citation = string (citation);
end
