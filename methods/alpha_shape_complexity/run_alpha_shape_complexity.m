% run alpha-shapes complexity code on SHREC collections
% james gardiner jdg@liverpool.ac.uk

close all
clear

% names of collections to process

all_collections = {'collection1','collection2_final','collection3'};

for i = 1:length(all_collections)
%     for i = 3
    
    clearvars -except i all_collections
    
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
        
        clearvars -except i j all_collections collection all_foldernames
        
        % for each subfolder generate a list of filenames
        
        foldername = all_foldernames{j};
        
        all_filenames = dir(fullfile('pre_processed',collection,...
            foldername,'*.mat'));
        
        for k = 1:length(all_filenames)
        % for k = 1 
            
            clearvars -except i j k all_collections collection all_foldernames foldername all_filenames
            close all
            
            tic
            
            % load in current files data
            
            filename = all_filenames(k).name;
            
            disp(['working on ' collection ' ' foldername ' ' filename])
            
            load(fullfile('pre_processed',collection,...
                foldername,filename))
            
            
            % calculate refinement coefficients and alpha radii to be used
            
            number_of_fits = 10;
            
            refinement_coefficients = exp(linspace(-0.8,5.7,number_of_fits));
            
            radii = refinement_coefficients.*scale; % scale is from the loaded data
            
            % set up variables for loop
            
            alpha_volumes = zeros(length(refinement_coefficients),1);
            plot_positions = [2 5 8];
            all_colors = hsv(100);
            color_to_plot = all_colors(randi(100),:);
            
            figure('color','w','position',[1 41 1280 607.3333])
            
            subplot_count = 1;
            
            % loop to run the alpha shapes code for each coefficient
            
            for m = 1:length(refinement_coefficients)
                
                if ismember(m,plot_positions) % if plotting result too
                    
                    subplot(1,length(plot_positions),subplot_count)
                    
                    [alpha_volumes(m),S] = alphavol([x y z],radii(m),0);
                    
                    h = trisurf(S.bnd,x,y,z);
                    set(h,'facecolor',color_to_plot,'edgecolor','none')
                    light
                    h.AmbientStrength = 0.3;
                    h.SpecularStrength = 0;
                    h.DiffuseStrength = 0.8;
                    axis equal
                    axis off
                    view(120,30)
                    
                    subplot_count = subplot_count +1;
                    
                else
                    
                    [alpha_volumes(m)] = alphavol([x y z],radii(m),0);
                    
                end
                
                % display progress
                
                if ~rem(m,round(number_of_fits/10))
                    
                    disp([num2str(m*100/number_of_fits) '%'])
                end
                
                
            end
            
            
            
            % create output folders if needed
            
            if ~isfolder(fullfile('figures','results',collection,foldername))
                
                mkdir(fullfile('figures','results',collection,foldername))
                
            end
            
            if ~isfolder(fullfile('results',collection,foldername))
                
                mkdir(fullfile('results',collection,foldername))
                
            end
            
            % write out figure and data
            
            drawnow
            
            print(fullfile('figures','results',collection,foldername,filename(1:end-4)),'-dpng')
            
            save(fullfile('results',collection,foldername,[filename(1:end-4) '.mat']),...
                'refinement_coefficients','alpha_volumes','scale','volume') % volume is from the loaded data
            
            toc
            
        end
        
    end
    
end