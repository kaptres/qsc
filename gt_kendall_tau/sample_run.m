% this outputs a table with hard-coded latex references and table alignments
all_scores;         % create the needed variables in the workspace

%% Collection 1
% additive perturbations
cellarr = col11;      % the first value in the table  (here, additively perturbed cubes)
cellarr2= col13;      % the second value in the table (here, additively perturbed spheres)

% subtractive perturbations (manually trim the rows of tables corresponding to the methods without any reported scores)
%cellarr = col12;      % the first value in the table  (here, subtractively perturbed cubes)
%cellarr2= col14;      % the second value in the table (here, subtractively perturbed spheres)
create_col1_table; 


%% Collection 2
%create_col2_table;


%% Collection 3
%create_col3_table;
