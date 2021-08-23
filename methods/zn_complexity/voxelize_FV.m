% author: mferhata
function vol = voxelize_FV (FV)
    % for more efficient use of memory
    Y_range = max(FV.vertices(:,1)) - min(FV.vertices(:,1));
    X_range = max(FV.vertices(:,2)) - min(FV.vertices(:,2));
    Z_range = max(FV.vertices(:,3)) - min(FV.vertices(:,3));
    target_volume = 300^3;
    const = (target_volume / (Y_range * X_range * Z_range)) ^ (1/3);

    vol = polygon2voxel (FV, floor (const * [X_range Y_range Z_range]), 'auto');
    vol = imfill (vol, 'holes');
end
