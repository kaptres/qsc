First compile buildScreenedPoissonSystemMex.c

mex -largeArrayDims buildScreenedPoissonSystemMex.c

compute_complexity function takes two parameters as:
dirPath : path of the directory containing 2D views of a 3D shape
shapeName : name of input 3D shape 
and returns complexity value of the input shape

Example run:
compute_complexity('c+', '000')
