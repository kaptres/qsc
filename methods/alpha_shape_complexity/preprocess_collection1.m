% Preprocess collection 1 data for SHREC 2021
% james gardiner jdg@liverpool.ac.uk

close all
clear

foldername = {'s+','s-','c+','c-'};

for j = 1:4
    
    clearvars -except foldername j
    
    all_filenames = dir(fullfile('raw_data','collection1_MATLAB',foldername{j},'*.mat'));
    
    
    for i = 1:length(all_filenames)
        % for i = 1
        
        close all
        clearvars -except all_filenames foldername i j
        
        % load in data
        
        filename = all_filenames(i).name;
        load(fullfile('raw_data','collection1_MATLAB',foldername{j},filename))
        
        disp([foldername{j} ' ' filename])
        
        
        % convert from indices to x, y, z coordinates
        
        [x_raw,y_raw,z_raw] = ind2sub(size(S),find(S));
        
        
        % calculate volume (i.e. in this case just number of elements)
        
        volume = length(x_raw);
        
        
        % downsample data
        
        downsample_level = 1e5;
        index = randsample(length(x_raw),downsample_level);
        x = x_raw(index);
        y = y_raw(index);
        z = z_raw(index);
        
        
        % calculate distnace to nearest neighbours to allow object to be scaled
        
        K = 100; % number of neighbours to calculate to.
        
        scale = nearest_neighbour_2([x y z],K,1e4);
        
        
        % plot to check
        
        figure('color','w','Position',[1 41 1280 6.073333333333333e+02])
        
        subplot(1,2,1)
        pcshow([x_raw y_raw z_raw])
        colormap('jet')
        axis equal
        axis off
        title('original point cloud')
        
        subplot(1,2,2)
        pcshow([x y z])
        colormap('jet')
        axis equal
        axis off
        title('downsampled point cloud')
        
        % write out figure and data
        
        drawnow
        
        print(fullfile('figures','pre_processed','collection1',foldername{j},filename(1:end-4)),'-dpng')
        
        % write out data
        
        save(fullfile('pre_processed','collection1',foldername{j},filename),'x','y','z','volume','scale')
        
    end
    
end