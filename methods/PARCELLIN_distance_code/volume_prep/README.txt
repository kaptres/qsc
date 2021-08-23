First compile polygon2voxel_double.c
mex -compatibleArrayDims polygon2voxel_double.c

In order to produce the complexity values presented in SHREC Track,
each 3D mesh model is voxelized to have around 300000-400000 voxels.
In this directory, the codes used for voxelization are presented.

The code used for voxelization is as follows:
[volume, numVoxels, vertex2Voxel] = voxelizeMesh(vertices, faces, scaleFactor, numVoxelsMin, numVoxelsMax)

numVoxelsMin and numVoxelsMax are selected as 300000 and 400000, respectively.
scaleFactor is selected as follows.

0.3 for c+, c-, s+, s-
2 for collection2_1
100 for collection2_2
200 for mesh segmentation benchmark

Example run:
fv = OBJ2FV(read_wobj('collection2_1_00.obj'));
[volume, numVoxels, vertex2Voxel] = voxelizeMesh(fv.vertices, fv.faces, 2, 300000, 400000);