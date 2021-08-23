% collate results for SHREC 2021
% james gardiner jdg@liverpool.ac.uk

close all
clear

%

line_count = 0;

% names of collections to process

all_collections = {'collection1','collection2_final','collection3'};

% for i = 1:length(all_collections)
for i = 1:3
    
    clearvars -except i all_collections all_data line_count
    
    % link collections to names of subfolders
    
    collection = all_collections{i};
    
    if strcmp(collection,'collection1')
        
        all_foldernames = {'s+','s-','c+','c-'};
        
    elseif strcmp(collection,'collection2_final')
        
        all_foldernames = {'collection2_family1','collection2_family2'};
        
    elseif strcmp(collection,'collection3')
        
        all_foldernames = {'off'};
        
    else
        
        disp(collection)
        warning('collection name not recognised')
        return
        
    end
    
    for j = 1:length(all_foldernames)
        % for j = 4
        
        clearvars -except i j all_collections collection all_foldernames all_data line_count
        
        % for each subfolder generate a list of filenames
        
        foldername = all_foldernames{j};
        
        all_filenames = dir(fullfile('pre_processed',collection,...
            foldername,'*.mat'));
        
        for k = 1:length(all_filenames)
            % for k = 1:3
            
            clearvars -except i j k all_collections collection all_foldernames foldername all_filenames all_data line_count
            close all
            
            
            % load in current files data
            
            filename = all_filenames(k).name;
            
            disp(['working on ' collection ' ' foldername ' ' filename])
            
            load(fullfile('results',collection,...
                foldername,filename))
            
            line_count = line_count + 1;
            
            all_data{line_count,1} = collection;
            all_data{line_count,2} = foldername;
            all_data{line_count,3} = filename(1:end-4);
            all_data{line_count,4} = scale;
            all_data{line_count,5} = volume;
            all_data{line_count,6} = alpha_volumes';
            
        end
        
    end
    
end

save(fullfile('results','collated_data'),'all_data','refinement_coefficients')
