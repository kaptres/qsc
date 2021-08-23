% run PCA on collated data
% james gardiner jdg@liverpool.ac.uk

close all
clear

% load in collated data

load(fullfile('results','collated_data'))


% chose which collection to run analysis on

run_all_collections = 0;

if run_all_collections
    
    collection_to_run = 'all';
    
else
    
    collection_to_run = 'collection1';
    % collection_to_run = 'collection2_final';
    % collection_to_run = 'collection3';
    
end

% collect data for selected collection

count = 0;

for i = 1:length(all_data)
    
    if run_all_collections
        
        count = count + 1;
        
        pca_inputs(count,:) = all_data{i,6}./all_data{i,5};
        
        pca_folders{count} = all_data{i,2};
        
        pca_models{count} = all_data{i,3};
        
    else
        
        if strcmp(all_data{i,1},collection_to_run)
            
            count = count + 1;
            
            pca_inputs(count,:) = all_data{i,6}./all_data{i,5};
            
            pca_folders{count} = all_data{i,2};
            
            pca_models{count} = all_data{i,3};
            
        end
        
    end
    
end

% sort collection3 to numerical order of models

if strcmp(collection_to_run,'collection3')
    
    for j = 1:length(pca_models)
        
        original_model_order(j) = str2num(pca_models{j});
        
    end
    
    original_pca_inputs = pca_inputs;
    original_pca_folders = pca_folders;
    original_pca_models = pca_models;
    
    [new_model_order,model_order_index] = sort(original_model_order);
    
    for j = 1:length(pca_models)
        
        pca_inputs(j,:) = original_pca_inputs(model_order_index(j),:);
        pca_folders{j} = original_pca_folders{model_order_index(j)};
        pca_models{j} = original_pca_models{model_order_index(j)};
        
    end
    
    
end



% run PCA

pca_inputs = zscore(pca_inputs); % to make equivalent to correlation matrix PCA
[coeff,score,latent,tsquared,explained] = pca(pca_inputs,'VariableWeights','variance');


for i = 1:length(refinement_coefficients)
    vbls{i} = num2str(refinement_coefficients(i));
end


unique_folders = unique(pca_folders);


% produce PCA biplot

% subplot(2,1,2)

figure('color','w')

for i = 1:length(unique_folders)
    
    for j = 1:length(pca_folders)
        
        index(j) = strcmp(pca_folders{j},unique_folders{i});
        
    end
    
    plot(score(index,1),score(index,2),'o')
    hold on
    
end

h = biplot(10*coeff(:,1:2),'varlabels',vbls);
set(h,'color','k')

h = legend(unique_folders);
set(h,'interpreter','none')

box off
xlabel(['PC1 (' num2str(explained(1),4) '%)'])
ylabel(['PC2 (' num2str(explained(2),4) '%)'])


% save PCA figure

drawnow

print(fullfile('figures','results',[collection_to_run '_PCA']),'-dpng')

% write out data in format needed for SHREC


if ~run_all_collections
    
    
    fid = fopen(fullfile('results',[collection_to_run '_output.txt']),'w');
    
    for i = 1:length(unique_folders)
        
        fprintf(fid,'%s\n',unique_folders{i});
        
        index = strcmp(pca_folders,unique_folders{i});
        
        for j = find(index)
            
            fprintf(fid,'(%f %f)',score(j,1),score(j,2));
            
        end
        
        fprintf(fid,'\n');
        
    end
    
    fclose(fid);
    
    
    % model check
    % same output as above expect model name instead of score just to ensure results are in proper order
    
    
    fid = fopen(fullfile('results',[collection_to_run '_output_model_check.txt']),'w');
    
    for i = 1:length(unique_folders)
        
        fprintf(fid,'%s\n',unique_folders{i});
        
        index = strcmp(pca_folders,unique_folders{i});
        
        for j = find(index)
            
            fprintf(fid,'(%s)',pca_models{j});
            
        end
        
        fprintf(fid,'\n');
        
    end
    
    fclose(fid);
    
end

