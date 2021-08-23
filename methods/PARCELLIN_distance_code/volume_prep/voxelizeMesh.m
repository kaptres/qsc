function [volume, numVoxels, vertex2Voxel] = voxelizeMesh(vertices, faces, scaleFactor, numVoxelsMin, numVoxelsMax)

fv.vertices = vertices;
fv.faces = faces;

% Process
sfactor = scaleFactor;
factStep = sfactor*0.05;
[volume,numVoxels,vertex2Voxel] = voxelize_factor(fv,sfactor);
while numVoxels > numVoxelsMax
	sfactor = sfactor-factStep;
	[volume,numVoxels,vertex2Voxel] = voxelize_factor(fv,sfactor);
	numVoxels
end
while numVoxels < numVoxelsMin
	sfactor = sfactor+factStep;
	[volume,numVoxels,vertex2Voxel] = voxelize_factor(fv,sfactor);
	numVoxels
end
numVoxels
end

function [volume,numVoxels,vertex2VoxelMap] = voxelize_factor(fv,scaleFactor)

% rescale integer grid resolution
fv.vertices = scaleFactor*fv.vertices;

% all vertices to have positive coordinates
fv.vertices(:,1) = fv.vertices(:,1)-min(fv.vertices(:,1))+1;
fv.vertices(:,2) = fv.vertices(:,2)-min(fv.vertices(:,2))+1;
fv.vertices(:,3) = fv.vertices(:,3)-min(fv.vertices(:,3))+1;

% determine ranges within individual axes
xrng = range(fv.vertices(:,1));
yrng = range(fv.vertices(:,2));
zrng = range(fv.vertices(:,3));

[volume] = polygon2voxel(fv,ceil([yrng xrng zrng])+2,'none');
volume = padarray(volume,[1 1 1]);
volume = imfill(volume,'holes');

numVoxels = nnz(volume);

vertex2VoxelMap = mapVertices2Voxels(fv.vertices,volume);
end

function vertex2VoxelMap = mapVertices2Voxels(vertices,volume)
% find correspondences bw mesh vertices and voxels
vertex2VoxelMap = round(vertices)+1;
corr = volume(sub2ind(size(volume),vertex2VoxelMap(:,2),vertex2VoxelMap(:,1),vertex2VoxelMap(:,3)));
if min(corr) == 0
	error('Nearest voxel is not included by volume');
end
end