% Preprocess collection 3 data for SHREC 2021
% james gardiner jdg@liverpool.ac.uk

close all
clear

foldername = 'off';
all_filenames = dir(fullfile('raw_data','collection3','MeshsegBenchmark-1.0','data',...
    foldername,'*.off'));

for i = 1:length(all_filenames)
% for i = 1
    
    tic
    
    close all
    clearvars -except all_filenames foldername i
    
    % load in data
    
    filename = all_filenames(i).name
    
    % load in off
    
    [str.vertices,str.faces] = read_off(fullfile('raw_data','collection3','MeshsegBenchmark-1.0','data',...
        foldername,filename));    
    
    
    % calculate volume ()
    
    [volume,~] = stlVolume(str.vertices,str.faces);
    
    if i == 261 % fix problem with the volume calculation for model 351.0ff
        
        volume = 0.229490;
        
    end
    
    
    % fill out internal volume of obj surface with points
    
    point_cloud_original = str.vertices';
    min_x = min(point_cloud_original(:,1));
    max_x = max(point_cloud_original(:,1));
    min_y = min(point_cloud_original(:,2));
    max_y = max(point_cloud_original(:,2));
    min_z = min(point_cloud_original(:,3));
    max_z = max(point_cloud_original(:,3));
    
    point_cloud = zeros(1.1e5,3); % empty matrix to fill with points
    
    while sum(point_cloud(:,1) ~= 0) < 1e5
        
        new_points(:,1) = min_x + (max_x - min_x).*rand(1e3,1);
        new_points(:,2) = min_y + (max_y - min_y).*rand(1e3,1);
        new_points(:,3) = min_z + (max_z - min_z).*rand(1e3,1);
        
        inside = in_polyhedron(str,new_points); % check which points are inside original mesh
        
        last_row = find(point_cloud(:,1) ~= 0,1,'last');
        
        if isempty(last_row) % first loop through find doesn't find anything
            
            last_row = 0;
            
        end
        
        all_rows = last_row+1:last_row+sum(inside);
        
        point_cloud(all_rows,:) = new_points(inside,:); % add new points to point cloud
        
    end
    
    point_cloud = point_cloud(1:1e5,:);  % cut down to 1e5 points
    
    x = point_cloud(:,1);
    y = point_cloud(:,2);
    z = point_cloud(:,3);
    
    
    % calculate distnace to nearest neighbours to allow object to be scaled
    
    K = 100; % number of neighbours to calculate to.
    
    scale = nearest_neighbour_2([x y z],K,1e4);
    
    
    % plot to check
    
    
    figure('color','w','Position',[1 41 1280 6.073333333333333e+02])
    
    subplot(1,3,1)
    h = trisurf(str.faces',str.vertices(1,:)',str.vertices(2,:)',str.vertices(3,:)');
    set(h,'facecolor','k','edgecolor','w')
    axis equal
    axis off
    h = title('original mesh');
    set(h,'color','w')
    
    subplot(1,3,2)
    pcshow(point_cloud_original)
    colormap('jet')
    axis equal
    axis off
    title('original point cloud')
    
    subplot(1,3,3)
    pcshow(point_cloud)
    colormap('jet')
    axis equal
    axis off
    title('filled point cloud')
    
    
    % write out figure and data
    
    drawnow
    
    print(fullfile('figures','pre_processed','collection3','off',filename(1:end-4)),'-dpng')
    
    save(fullfile('pre_processed','collection3','off',[filename(1:end-4) '.mat']),'x','y','z','volume','scale')
    
    toc
    
end

