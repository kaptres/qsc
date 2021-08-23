function [mean_distance_nearest_neighbour] = nearest_neighbour_2(point_cloud,K,N)

% james gardiner jdg@liverpool.ac.uk

% calculate distance to nearest neighbours as an estimate of point cloud
% density. K is the number of neighbours to use. N is number of points to
% perform calculation for.

ptCloud = pointCloud(point_cloud);

dists = zeros(N,K);
count = 0;

for i = randi(length(point_cloud(:,1)),[1 N])
    
    count = count+1;
    
    [~,dists(count,:)] = findNearestNeighbors(ptCloud,ptCloud.Location(i,:),K);
    
end

mean_distance_nearest_neighbour = mean(mean(dists));



