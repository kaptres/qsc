% Preprocess collection 2 data for SHREC 2021
% james gardiner jdg@liverpool.ac.uk

close all
clear


% choose folder to run on

foldername = 'collection2_family1';
% foldername = 'collection2_family2';


all_filenames = dir(fullfile('raw_data','collection2_final',...
    foldername,'*.obj'));

for i = 1:length(all_filenames)
    % for i = 1
    
    tic
    
    close all
    clearvars -except all_filenames foldername i
    
    % load in data
    
    filename = all_filenames(i).name
    
    % load in obj
    
    obj = readObj(fullfile('raw_data','collection2_final',...
        foldername,filename));
    
    str.vertices = obj.v;
    str.faces = obj.f.v;
    
    % calculate volume ()
    
    [volume,~] = stlVolume(str.vertices',str.faces');
    
    
    % fill out internal volume of obj surface with points
    
    point_cloud_original = str.vertices;
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
    h = trisurf(str.faces,str.vertices(:,1),str.vertices(:,2),str.vertices(:,3));
    set(h,'facecolor','k','edgecolor','w')
    axis equal
    axis off
    h = title('original mesh');
    set(h,'color','w')
    
    subplot(1,3,2)
    h = pcshow(point_cloud_original);
    h2 = h.Children;
    set(h2,'SizeData',20)
    colormap('jet')
    axis equal
    axis off
    title('original point cloud')
    
    subplot(1,3,3)
    h = pcshow(point_cloud);
    h2 = h.Children;
    set(h2,'SizeData',20)
    colormap('jet')
    axis equal
    axis off
    title('filled point cloud')
    
    
    % write out figure and data
    
    drawnow
    
    print(fullfile('figures','pre_processed','collection2_final',foldername,filename(1:end-4)),'-dpng')
    
    save(fullfile('pre_processed','collection2_final',foldername,[filename(1:end-4) '.mat']),'x','y','z','volume','scale')
    
    toc
    
end

